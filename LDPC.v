// LDPC.v -- Lab02 v6_noB2 architecture.
// v1 base: 16 CNU lanes (one full layer per cycle), rotated storage in
// bank A/B, touch-based (first/mid/last) bank control shared by both
// modes, 4-stage c2v FIFO, rotate-by-1 input load / output.
// v2 added (arch_notes.md 4.1, unchanged numerics / latency):
//   O1: edge arithmetic in sign-magnitude form (op/carry-in adds instead
//       of negate-then-subtract; merged clip+abs bit-trick; CNU ports
//       take/return magnitude+sign directly, no abs6/make_r).
//   O2: the layer-dependent read-side edge table (col_sel/touch_first,
//       replicated x4 for fanout) and every per-column bank-writeback
//       select are decoded one cycle ahead into registers instead of
//       live-decoding the `layer` counter each cycle.
// v3 added (arch_notes.md 4.2, latency 4540 -> 4340):
//   Core: the last input word (Lch127, needed by lane 12 / column 7 /
//       edge 6 at iteration-0 layer-0) is folded into LOAD's last cycle
//       instead of costing its own RUN cycle, via a capture register
//       (cap_desc) plus a short fast path combined with in_data.
//   S1/S2: norm5 as a 32-entry table; qsel registered like touch_first.
// v4 adds (arch_notes.md 4.3, unchanged numerics / latency):
//   S6: B drops its "hold" leg for every column except 1/2 (the only two
//       that hit an absent layer while their B value is still live) --
//       those columns' B write is now an unconditional 2/3-way select
//       among its real sources, never gated by run_now.
//   S7: slot 1 (edge 0)'s B side is hardwired to B[1] -- column 0's B is
//       never actually selected there (column 0's only slot-1 touch is
//       its first touch, always reading A), so col_sel never needs to
//       reach B at all for this edge.
//   S8: A drops its hold too. Every column now does *something* every
//       cycle: LOAD/OUTPUT/STOP/IDLE all rotate by 1 (reusing the same
//       leg instead of a dedicated self-feedback mux), LOAD's injection
//       point moves accordingly (inj[c] = (-tf[c]-2) mod 16), and each
//       column's absent-layer action is either a free reuse of an
//       existing delta-0 leg (columns 0/1) or an explicit rot1 (columns
//       2/3) depending on which is cheaper. The capture step (for
//       lane 12's edge 6 at the overlap cycle) now reads/forces lane 13
//       instead of 12, since column 1-6 keep rotating through the whole
//       LOAD window. See arch5.py / search3.py for the derivation.
// v5 (arch_notes.md 4.4, unchanged numerics / latency):
//   overlap_cycle no longer reads in_data_valid -- see the comment at
//   its definition below. This keeps in_data_valid's 0.5T input delay
//   off the run_now/commit path that fans out to every A/B write select,
//   which was the first slack-0 path found ahead of a T sweep.
//   (A) lane 12's fast path shares its two nx comparisons (8 -> 2).
// v6 adds (arch_notes.md 4.6, unchanged numerics / latency):
//   v6-1 (key CNU): each edge's exact two's-complement magnitude (which
//       needed a per-edge +1 incrementer to turn one's-complement into
//       two's-complement) is replaced by a 6-bit key K = {one's-complement
//       magnitude, sign} that CNU_lane's tree compares directly. The
//       deferred +1 correction and the (3m+2)>>2 normalize are folded into
//       one 64-entry table (nk64) applied only at the two root outputs
//       instead of once per edge -- 112 per-edge incrementers become 32
//       table lookups (16 lanes x 2 outputs).
//   v6-2 (B1): the overlap cycle's fast path no longer waits for
//       norm(|in_data|) to gate its two comparisons against cap_desc's
//       cn1/cn2 -- since norm is monotonic, nx<c is equivalent to
//       |in_data|<thr(c) for a precomputed threshold table (thr24),
//       which only depends on cap_desc (stable well before in_data
//       arrives). in_data is compared directly against +-thr instead of
//       through abs+norm first. Also drops the now-unneeded crmag mux
//       (edge 6's capture-forced key can never win the tie, so its fast
//       path value is simply cn1).
// v6_noB2 keeps v6-1 and v6-2 but removes v6-3 (B2): CNU_lane's
// final merge now selects the raw min1/min2 keys first and applies nk64
// only to those two root outputs. This is numerically equivalent and
// trades root-path timing for lower CNU area.
// See arch_notes.md section 3 (base architecture) and 4.2/4.3/4.4/4.6
// (v3/v4/v5/v6 derivation) for the numeric derivation of every constant
// below.

module CNU_lane (
    input      [5:0] key0,
    input      [5:0] key1,
    input      [5:0] key2,
    input      [5:0] key3,
    input      [5:0] key4,
    input      [5:0] key5,
    input      [5:0] key6,
    input      [6:0] qneg,

    output reg [19:0] c2v_new,
    output reg [34:0] rmag_flat,
    output reg [6:0]  rneg
);

    reg [14:0] leaf [0:6];
    reg [14:0] pair01;
    reg [14:0] pair23;
    reg [14:0] pair45;
    reg [14:0] group03;
    reg [14:0] group46;
    reg [14:0] root_pair;

    reg [4:0] norm_min1;
    reg [4:0] norm_min2;
    reg [2:0] root_idx;
    reg       sign_parity;

    // key = {one's-complement magnitude (5b, saturated), sign}. Ties
    // always resolve to the lower edge index (see merge_top2 below), so a
    // forced key of 63 (the maximum possible 6-bit value) can never win --
    // that is what lets the capture path (LANE13_SLOT7) force-lose edge 6.
    function [14:0] merge_top2;
        input [14:0] left_value;
        input [14:0] right_value;
        reg [5:0] next_min2;
        begin
            // The tree always places lower edge indices in the left subtree.
            // Selecting left on a tie therefore implements the required tie rule.
            if (left_value[14:9] <= right_value[14:9]) begin
                if (left_value[5:0] < right_value[14:9])
                    next_min2 = left_value[5:0];
                else
                    next_min2 = right_value[14:9];
                merge_top2 = {left_value[14:6], next_min2};
            end else begin
                if (right_value[5:0] < left_value[14:9])
                    next_min2 = right_value[5:0];
                else
                    next_min2 = left_value[14:9];
                merge_top2 = {right_value[14:6], next_min2};
            end
        end
    endfunction

    // v6-1: nk64 folds two steps into one 64-entry table on the 6-bit key:
    // (1) the one's-complement -> two's-complement +1 correction that used
    // to run once per edge (mp' = mp + sign, saturate at 31), and (2) S1's
    // (3m+2)>>2 normalize. Equivalent to norm5(min((K>>1)+(K&1),31)).
    // Applied only at the two root outputs (32 lookups total) instead of
    // once per edge (112 incrementers removed).
    function [4:0] nk64;
        input [5:0] k;
        begin
            case (k)
                6'd0:  nk64 = 5'd0;
                6'd1:  nk64 = 5'd1;
                6'd2:  nk64 = 5'd1;
                6'd3:  nk64 = 5'd2;
                6'd4:  nk64 = 5'd2;
                6'd5:  nk64 = 5'd2;
                6'd6:  nk64 = 5'd2;
                6'd7:  nk64 = 5'd3;
                6'd8:  nk64 = 5'd3;
                6'd9:  nk64 = 5'd4;
                6'd10: nk64 = 5'd4;
                6'd11: nk64 = 5'd5;
                6'd12: nk64 = 5'd5;
                6'd13: nk64 = 5'd5;
                6'd14: nk64 = 5'd5;
                6'd15: nk64 = 5'd6;
                6'd16: nk64 = 5'd6;
                6'd17: nk64 = 5'd7;
                6'd18: nk64 = 5'd7;
                6'd19: nk64 = 5'd8;
                6'd20: nk64 = 5'd8;
                6'd21: nk64 = 5'd8;
                6'd22: nk64 = 5'd8;
                6'd23: nk64 = 5'd9;
                6'd24: nk64 = 5'd9;
                6'd25: nk64 = 5'd10;
                6'd26: nk64 = 5'd10;
                6'd27: nk64 = 5'd11;
                6'd28: nk64 = 5'd11;
                6'd29: nk64 = 5'd11;
                6'd30: nk64 = 5'd11;
                6'd31: nk64 = 5'd12;
                6'd32: nk64 = 5'd12;
                6'd33: nk64 = 5'd13;
                6'd34: nk64 = 5'd13;
                6'd35: nk64 = 5'd14;
                6'd36: nk64 = 5'd14;
                6'd37: nk64 = 5'd14;
                6'd38: nk64 = 5'd14;
                6'd39: nk64 = 5'd15;
                6'd40: nk64 = 5'd15;
                6'd41: nk64 = 5'd16;
                6'd42: nk64 = 5'd16;
                6'd43: nk64 = 5'd17;
                6'd44: nk64 = 5'd17;
                6'd45: nk64 = 5'd17;
                6'd46: nk64 = 5'd17;
                6'd47: nk64 = 5'd18;
                6'd48: nk64 = 5'd18;
                6'd49: nk64 = 5'd19;
                6'd50: nk64 = 5'd19;
                6'd51: nk64 = 5'd20;
                6'd52: nk64 = 5'd20;
                6'd53: nk64 = 5'd20;
                6'd54: nk64 = 5'd20;
                6'd55: nk64 = 5'd21;
                6'd56: nk64 = 5'd21;
                6'd57: nk64 = 5'd22;
                6'd58: nk64 = 5'd22;
                6'd59: nk64 = 5'd23;
                6'd60: nk64 = 5'd23;
                default: nk64 = 5'd23; // k == 61, 62, 63
            endcase
        end
    endfunction

    always @(*) begin
        leaf[0] = {key0, 3'd0, 6'd63};
        leaf[1] = {key1, 3'd1, 6'd63};
        leaf[2] = {key2, 3'd2, 6'd63};
        leaf[3] = {key3, 3'd3, 6'd63};
        leaf[4] = {key4, 3'd4, 6'd63};
        leaf[5] = {key5, 3'd5, 6'd63};
        leaf[6] = {key6, 3'd6, 6'd63};

        pair01  = merge_top2(leaf[0], leaf[1]);
        pair23  = merge_top2(leaf[2], leaf[3]);
        pair45  = merge_top2(leaf[4], leaf[5]);
        group03 = merge_top2(pair01, pair23);
        group46 = merge_top2(pair45, leaf[6]);

        // v6_noB2: merge raw keys before normalizing the two root minima.
        root_pair = merge_top2(group03, group46);
        root_idx  = root_pair[8:6];
        norm_min1 = nk64(root_pair[14:9]);
        norm_min2 = nk64(root_pair[5:0]);

        sign_parity = ^qneg;
        rneg = qneg ^ {7{sign_parity}};

        c2v_new = {norm_min1, norm_min2, root_idx, rneg};

        rmag_flat[4:0]   = (root_idx == 3'd0) ? norm_min2 : norm_min1;
        rmag_flat[9:5]   = (root_idx == 3'd1) ? norm_min2 : norm_min1;
        rmag_flat[14:10] = (root_idx == 3'd2) ? norm_min2 : norm_min1;
        rmag_flat[19:15] = (root_idx == 3'd3) ? norm_min2 : norm_min1;
        rmag_flat[24:20] = (root_idx == 3'd4) ? norm_min2 : norm_min1;
        rmag_flat[29:25] = (root_idx == 3'd5) ? norm_min2 : norm_min1;
        rmag_flat[34:30] = (root_idx == 3'd6) ? norm_min2 : norm_min1;
    end

endmodule


module LDPC (
    input              clk,
    input              rst_n,
    input              in_mode_valid,
    input              in_mode,
    input              in_data_valid,
    input      [5:0]   in_data,
    output reg         out_valid,
    output reg [7:0]   out_data,
    output reg         out_warn
);

    localparam [1:0] S_IDLE = 2'd0;
    localparam [1:0] S_LOAD = 2'd1;
    localparam [1:0] S_RUN  = 2'd2;
    localparam [1:0] S_OUT  = 2'd3;

    reg [1:0] state;
    reg [1:0] next_state;
    reg       mode_reg;

    reg [6:0] io_cnt;
    reg [1:0] layer;
    reg [3:0] iter;

    reg signed [7:0] A [0:7][0:15];
    reg signed [7:0] B [0:7][0:15];
    reg        [19:0] c2v_fifo [0:3][0:15];

    integer rr;

    //============================================================
    // Shared control
    //============================================================

    wire check_cycle = (state == S_RUN) && (layer == 2'd0) && (iter != 4'd0);
    wire at_limit     = (iter == 4'd8);
    wire syndrome_nonzero;
    wire stop_now = check_cycle && (!syndrome_nonzero || at_limit);

    // The overlap cycle is LOAD's very last cycle (io_cnt==127); it
    // commits iteration-0 layer-0 for every lane exactly like a RUN cycle.
    // v5: dropped the in_data_valid term. PATTERN guarantees in_data_valid
    // is high for exactly 128 consecutive cycles and io_cnt only advances
    // while it is high, so "LOAD and io_cnt==127" already means this is
    // that last input cycle. Not reading in_data_valid here keeps its
    // 0.5T input delay off the overlap_cycle -> run_now -> commit path
    // that fans out to every A/B write select (this path was slack-0 at
    // T=10 and would be the first to fail as T sweeps down).
    wire overlap_cycle = (state == S_LOAD) && (io_cnt == 7'd127);
    wire run_now        = (state == S_RUN) || overlap_cycle;
    wire commit          = run_now && !stop_now;
    wire entering_run    = (next_state == S_RUN) && (state != S_RUN);

    always @(*) begin
        next_state = state;
        case (state)
            S_IDLE: if (in_mode_valid && !out_valid) next_state = S_LOAD;
            S_LOAD: if (in_data_valid && io_cnt == 7'd127) next_state = S_RUN;
            S_RUN:  if (stop_now) next_state = S_OUT;
            S_OUT:  if (io_cnt == 7'd127) next_state = S_IDLE;
            default: next_state = S_IDLE;
        endcase
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) state <= S_IDLE;
        else        state <= next_state;
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n)
            mode_reg <= 1'b0;
        else if (in_mode_valid && (state == S_IDLE) && !out_valid)
            mode_reg <= in_mode;
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            io_cnt <= 7'd0;
        end else if (in_mode_valid && (state == S_IDLE) && !out_valid) begin
            io_cnt <= 7'd0;
        end else if ((state == S_LOAD) && in_data_valid) begin
            io_cnt <= io_cnt + 7'd1;
        end else if (stop_now) begin
            io_cnt <= 7'd1;
        end else if (state == S_OUT) begin
            io_cnt <= (io_cnt == 7'd127) ? 7'd0 : io_cnt + 7'd1;
        end
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            layer <= 2'd0;
            iter  <= 4'd0;
        end else if (in_mode_valid && (state == S_IDLE) && !out_valid) begin
            layer <= 2'd0;
            iter  <= 4'd0;
        end else if (entering_run) begin
            // Layer 0 of iteration 1 was already committed during the
            // overlap cycle (still state S_LOAD); RUN itself starts at
            // layer 1.
            layer <= 2'd1;
        end else if ((state == S_RUN) && !stop_now) begin
            layer <= (layer == 2'd3) ? 2'd0 : layer + 2'd1;
            if (layer == 2'd3) iter <= iter + 4'd1;
        end
    end

    //============================================================
    // O2: next_layer_comb predicts, one cycle ahead, the layer value
    // that every layer-dependent select below should be decoded for.
    //   - io_cnt==126 (one cycle before the overlap cycle) -> layer 0
    //   - the overlap cycle itself (io_cnt==127)            -> layer 1
    //   - normal RUN advance                                 -> layer+1
    //   - otherwise (prediction unused)                      -> hold
    //============================================================

    wire pre_layer0 = (state == S_LOAD) && in_data_valid && (io_cnt == 7'd126);
    wire [1:0] next_layer_comb =
        pre_layer0                        ? 2'd0 :
        overlap_cycle                      ? 2'd1 :
        ((state == S_RUN) && !stop_now)   ? ((layer == 2'd3) ? 2'd0 : layer + 2'd1) :
                                             layer;

    // Shared one-hot mirror of "current layer", registered one cycle
    // ahead so the per-column bank-writeback blocks further below never
    // decode `layer` live; they just read nl_r0..nl_r3.
    reg nl_r0, nl_r1, nl_r2, nl_r3;
    always @(posedge clk) begin
        nl_r0 <= (next_layer_comb == 2'd0);
        nl_r1 <= (next_layer_comb == 2'd1);
        nl_r2 <= (next_layer_comb == 2'd2);
        nl_r3 <= (next_layer_comb == 2'd3);
    end

    //============================================================
    // Per-layer edge table: edge e (0..6) is normally column e+1;
    // column 0 borrows edge (layer-1) at layer 1/2/3.
    // Replicated into 4 groups (4 lanes each) so col_sel/touch_first
    // never fan out past the group of lanes they drive.
    // S2: qsel (= !mode | touch_first) is registered the same way,
    // so the Q mux never reads `mode_reg` live either.
    //============================================================

    reg [2:0] col_sel_g0 [0:6]; reg touch_first_g0 [0:6]; reg qsel_g0 [0:6];
    reg [2:0] col_sel_g1 [0:6]; reg touch_first_g1 [0:6]; reg qsel_g1 [0:6];
    reg [2:0] col_sel_g2 [0:6]; reg touch_first_g2 [0:6]; reg qsel_g2 [0:6];
    reg [2:0] col_sel_g3 [0:6]; reg touch_first_g3 [0:6]; reg qsel_g3 [0:6];
    integer ti;

    always @(posedge clk) begin
        case (next_layer_comb)
            2'd0: begin
                col_sel_g0[0] <= 3'd1; col_sel_g0[1] <= 3'd2; col_sel_g0[2] <= 3'd3;
                col_sel_g0[3] <= 3'd4; col_sel_g0[4] <= 3'd5; col_sel_g0[5] <= 3'd6; col_sel_g0[6] <= 3'd7;
                col_sel_g1[0] <= 3'd1; col_sel_g1[1] <= 3'd2; col_sel_g1[2] <= 3'd3;
                col_sel_g1[3] <= 3'd4; col_sel_g1[4] <= 3'd5; col_sel_g1[5] <= 3'd6; col_sel_g1[6] <= 3'd7;
                col_sel_g2[0] <= 3'd1; col_sel_g2[1] <= 3'd2; col_sel_g2[2] <= 3'd3;
                col_sel_g2[3] <= 3'd4; col_sel_g2[4] <= 3'd5; col_sel_g2[5] <= 3'd6; col_sel_g2[6] <= 3'd7;
                col_sel_g3[0] <= 3'd1; col_sel_g3[1] <= 3'd2; col_sel_g3[2] <= 3'd3;
                col_sel_g3[3] <= 3'd4; col_sel_g3[4] <= 3'd5; col_sel_g3[5] <= 3'd6; col_sel_g3[6] <= 3'd7;
                for (ti = 0; ti < 7; ti = ti + 1) begin
                    touch_first_g0[ti] <= 1'b1; touch_first_g1[ti] <= 1'b1;
                    touch_first_g2[ti] <= 1'b1; touch_first_g3[ti] <= 1'b1;
                    qsel_g0[ti] <= 1'b1; qsel_g1[ti] <= 1'b1;
                    qsel_g2[ti] <= 1'b1; qsel_g3[ti] <= 1'b1;
                end
            end
            2'd1: begin
                col_sel_g0[0] <= 3'd0; col_sel_g0[1] <= 3'd2; col_sel_g0[2] <= 3'd3;
                col_sel_g0[3] <= 3'd4; col_sel_g0[4] <= 3'd5; col_sel_g0[5] <= 3'd6; col_sel_g0[6] <= 3'd7;
                col_sel_g1[0] <= 3'd0; col_sel_g1[1] <= 3'd2; col_sel_g1[2] <= 3'd3;
                col_sel_g1[3] <= 3'd4; col_sel_g1[4] <= 3'd5; col_sel_g1[5] <= 3'd6; col_sel_g1[6] <= 3'd7;
                col_sel_g2[0] <= 3'd0; col_sel_g2[1] <= 3'd2; col_sel_g2[2] <= 3'd3;
                col_sel_g2[3] <= 3'd4; col_sel_g2[4] <= 3'd5; col_sel_g2[5] <= 3'd6; col_sel_g2[6] <= 3'd7;
                col_sel_g3[0] <= 3'd0; col_sel_g3[1] <= 3'd2; col_sel_g3[2] <= 3'd3;
                col_sel_g3[3] <= 3'd4; col_sel_g3[4] <= 3'd5; col_sel_g3[5] <= 3'd6; col_sel_g3[6] <= 3'd7;
                touch_first_g0[0] <= 1'b1; touch_first_g1[0] <= 1'b1;
                touch_first_g2[0] <= 1'b1; touch_first_g3[0] <= 1'b1;
                qsel_g0[0] <= 1'b1; qsel_g1[0] <= 1'b1; qsel_g2[0] <= 1'b1; qsel_g3[0] <= 1'b1;
                for (ti = 1; ti < 7; ti = ti + 1) begin
                    touch_first_g0[ti] <= 1'b0; touch_first_g1[ti] <= 1'b0;
                    touch_first_g2[ti] <= 1'b0; touch_first_g3[ti] <= 1'b0;
                    qsel_g0[ti] <= !mode_reg; qsel_g1[ti] <= !mode_reg;
                    qsel_g2[ti] <= !mode_reg; qsel_g3[ti] <= !mode_reg;
                end
            end
            2'd2: begin
                col_sel_g0[0] <= 3'd1; col_sel_g0[1] <= 3'd0; col_sel_g0[2] <= 3'd3;
                col_sel_g0[3] <= 3'd4; col_sel_g0[4] <= 3'd5; col_sel_g0[5] <= 3'd6; col_sel_g0[6] <= 3'd7;
                col_sel_g1[0] <= 3'd1; col_sel_g1[1] <= 3'd0; col_sel_g1[2] <= 3'd3;
                col_sel_g1[3] <= 3'd4; col_sel_g1[4] <= 3'd5; col_sel_g1[5] <= 3'd6; col_sel_g1[6] <= 3'd7;
                col_sel_g2[0] <= 3'd1; col_sel_g2[1] <= 3'd0; col_sel_g2[2] <= 3'd3;
                col_sel_g2[3] <= 3'd4; col_sel_g2[4] <= 3'd5; col_sel_g2[5] <= 3'd6; col_sel_g2[6] <= 3'd7;
                col_sel_g3[0] <= 3'd1; col_sel_g3[1] <= 3'd0; col_sel_g3[2] <= 3'd3;
                col_sel_g3[3] <= 3'd4; col_sel_g3[4] <= 3'd5; col_sel_g3[5] <= 3'd6; col_sel_g3[6] <= 3'd7;
                for (ti = 0; ti < 7; ti = ti + 1) begin
                    touch_first_g0[ti] <= 1'b0; touch_first_g1[ti] <= 1'b0;
                    touch_first_g2[ti] <= 1'b0; touch_first_g3[ti] <= 1'b0;
                    qsel_g0[ti] <= !mode_reg; qsel_g1[ti] <= !mode_reg;
                    qsel_g2[ti] <= !mode_reg; qsel_g3[ti] <= !mode_reg;
                end
            end
            default: begin // layer 3
                col_sel_g0[0] <= 3'd1; col_sel_g0[1] <= 3'd2; col_sel_g0[2] <= 3'd0;
                col_sel_g0[3] <= 3'd4; col_sel_g0[4] <= 3'd5; col_sel_g0[5] <= 3'd6; col_sel_g0[6] <= 3'd7;
                col_sel_g1[0] <= 3'd1; col_sel_g1[1] <= 3'd2; col_sel_g1[2] <= 3'd0;
                col_sel_g1[3] <= 3'd4; col_sel_g1[4] <= 3'd5; col_sel_g1[5] <= 3'd6; col_sel_g1[6] <= 3'd7;
                col_sel_g2[0] <= 3'd1; col_sel_g2[1] <= 3'd2; col_sel_g2[2] <= 3'd0;
                col_sel_g2[3] <= 3'd4; col_sel_g2[4] <= 3'd5; col_sel_g2[5] <= 3'd6; col_sel_g2[6] <= 3'd7;
                col_sel_g3[0] <= 3'd1; col_sel_g3[1] <= 3'd2; col_sel_g3[2] <= 3'd0;
                col_sel_g3[3] <= 3'd4; col_sel_g3[4] <= 3'd5; col_sel_g3[5] <= 3'd6; col_sel_g3[6] <= 3'd7;
                for (ti = 0; ti < 7; ti = ti + 1) begin
                    touch_first_g0[ti] <= 1'b0; touch_first_g1[ti] <= 1'b0;
                    touch_first_g2[ti] <= 1'b0; touch_first_g3[ti] <= 1'b0;
                    qsel_g0[ti] <= !mode_reg; qsel_g1[ti] <= !mode_reg;
                    qsel_g2[ti] <= !mode_reg; qsel_g3[ti] <= !mode_reg;
                end
            end
        endcase
    end

    //============================================================
    // Edge datapath (7 edges x 16 lanes).
    // O1: r_old subtraction done as a single add with conditional
    // invert + carry-in (op/op_s) instead of negate-then-subtract;
    // clip(|.|,31) computed directly in sign-magnitude form.
    // v3/v4: edge 6 (column 7) at lane 12 is special for x_base/newv
    // during the overlap cycle; at lane 13 it is special throughout
    // LOAD (forced) for the capture register below.
    //============================================================

    wire signed [7:0] Xv     [0:6][0:15];
    wire signed [7:0] Qv     [0:6][0:15];
    wire        [4:0] ro_mag [0:6][0:15];
    wire               ro_neg [0:6][0:15];
    wire signed [7:0] q_base [0:6][0:15];
    wire signed [7:0] x_base [0:6][0:15];
    // v6-1: q_key replaces q_mag -- one's-complement magnitude (no +1
    // correction) packed with the sign into a 6-bit key for CNU_lane's
    // tree. q_sgn is kept as a separate signal (same value as q_key's LSB
    // in the normal case) because the capture force below must decouple
    // them: key forces to 63 (always-loses) while q_sgn forces to 0.
    wire        [5:0] q_key  [0:6][0:15];
    wire               q_sgn  [0:6][0:15];

    genvar ge, gl;
    generate
        // Xv/Qv: 4 lane-groups, each fed from its own replicated
        // col_sel_gX / touch_first_gX / qsel_gX registers. S7: edge 0's
        // B side is hardwired to B[1] -- column 0's B is never selected
        // there (its only slot-1 touch is 'first', always reading A).
        for (ge = 0; ge < 7; ge = ge + 1) begin: XQ0
            for (gl = 0; gl < 4; gl = gl + 1) begin: LANE
                if (ge == 0) begin: EDGE0_S7
                    assign Xv[ge][gl] = touch_first_g0[ge] ? A[col_sel_g0[ge]][gl] : B[1][gl];
                    assign Qv[ge][gl] = qsel_g0[ge]        ? A[col_sel_g0[ge]][gl] : B[1][gl];
                end else begin: EDGE_NORMAL
                    assign Xv[ge][gl] = touch_first_g0[ge] ? A[col_sel_g0[ge]][gl] : B[col_sel_g0[ge]][gl];
                    assign Qv[ge][gl] = qsel_g0[ge]        ? A[col_sel_g0[ge]][gl] : B[col_sel_g0[ge]][gl];
                end
            end
        end
        for (ge = 0; ge < 7; ge = ge + 1) begin: XQ1
            for (gl = 4; gl < 8; gl = gl + 1) begin: LANE
                if (ge == 0) begin: EDGE0_S7
                    assign Xv[ge][gl] = touch_first_g1[ge] ? A[col_sel_g1[ge]][gl] : B[1][gl];
                    assign Qv[ge][gl] = qsel_g1[ge]        ? A[col_sel_g1[ge]][gl] : B[1][gl];
                end else begin: EDGE_NORMAL
                    assign Xv[ge][gl] = touch_first_g1[ge] ? A[col_sel_g1[ge]][gl] : B[col_sel_g1[ge]][gl];
                    assign Qv[ge][gl] = qsel_g1[ge]        ? A[col_sel_g1[ge]][gl] : B[col_sel_g1[ge]][gl];
                end
            end
        end
        for (ge = 0; ge < 7; ge = ge + 1) begin: XQ2
            for (gl = 8; gl < 12; gl = gl + 1) begin: LANE
                if (ge == 0) begin: EDGE0_S7
                    assign Xv[ge][gl] = touch_first_g2[ge] ? A[col_sel_g2[ge]][gl] : B[1][gl];
                    assign Qv[ge][gl] = qsel_g2[ge]        ? A[col_sel_g2[ge]][gl] : B[1][gl];
                end else begin: EDGE_NORMAL
                    assign Xv[ge][gl] = touch_first_g2[ge] ? A[col_sel_g2[ge]][gl] : B[col_sel_g2[ge]][gl];
                    assign Qv[ge][gl] = qsel_g2[ge]        ? A[col_sel_g2[ge]][gl] : B[col_sel_g2[ge]][gl];
                end
            end
        end
        for (ge = 0; ge < 7; ge = ge + 1) begin: XQ3
            for (gl = 12; gl < 16; gl = gl + 1) begin: LANE
                if (ge == 0) begin: EDGE0_S7
                    assign Xv[ge][gl] = touch_first_g3[ge] ? A[col_sel_g3[ge]][gl] : B[1][gl];
                    assign Qv[ge][gl] = qsel_g3[ge]        ? A[col_sel_g3[ge]][gl] : B[1][gl];
                end else begin: EDGE_NORMAL
                    assign Xv[ge][gl] = touch_first_g3[ge] ? A[col_sel_g3[ge]][gl] : B[col_sel_g3[ge]][gl];
                    assign Qv[ge][gl] = qsel_g3[ge]        ? A[col_sel_g3[ge]][gl] : B[col_sel_g3[ge]][gl];
                end
            end
        end

        for (ge = 0; ge < 7; ge = ge + 1) begin: EDGE
            for (gl = 0; gl < 16; gl = gl + 1) begin: LANE
                assign ro_mag[ge][gl] = (c2v_fifo[0][gl][9:7] == ge[2:0])
                                        ? c2v_fifo[0][gl][14:10]
                                        : c2v_fifo[0][gl][19:15];
                assign ro_neg[ge][gl] = c2v_fifo[0][gl][ge];

                // q_base = Q - r_old, via conditional-invert add.
                wire       op_s_w = ~ro_neg[ge][gl];
                wire [7:0] op_w   = {3'b0, ro_mag[ge][gl]} ^ {8{op_s_w}};
                assign q_base[ge][gl] = Qv[ge][gl] + op_w + {7'b0, op_s_w};

                if (ge == 6 && gl == 12) begin: LANE12_SLOT7
                    // Overlap cycle: Lch127 bypasses A_7[12] (garbage,
                    // since column 7's injection point never targets it);
                    // r_old is 0 at iteration 0 layer 0, so x_base = in_data.
                    wire signed [7:0] indata_ext = {{2{in_data[5]}}, in_data};
                    assign x_base[ge][gl] = overlap_cycle ? indata_ext
                                                            : (Xv[ge][gl] + op_w + {7'b0, op_s_w});

                    wire        sgn_w = q_base[ge][gl][7];
                    wire [6:0]  lo_w  = sgn_w ? ~q_base[ge][gl][6:0] : q_base[ge][gl][6:0];
                    wire        sat_w = (lo_w[6:5] != 2'b00);
                    wire [4:0]  mp_w  = sat_w ? 5'd31 : lo_w[4:0];
                    assign q_key[ge][gl] = {mp_w, sgn_w};
                    assign q_sgn[ge][gl] = sgn_w;
                end else if (ge == 6 && gl == 13) begin: LANE13_SLOT7
                    // Capture cycle input: forced to (mag 31, sign +) for
                    // the whole LOAD window so this edge never wins the
                    // CNU's min-search -- that is exactly what lets
                    // cap_desc (below) capture the other 6 real edges'
                    // stats. Column 1-6 keep rotating through all of LOAD
                    // now (S8), so by the time it matters (cycle 126) the
                    // stale, never-injected value has rotated to lane 13
                    // (one step ahead of lane 12, which is where the
                    // overlap cycle itself will look).
                    assign x_base[ge][gl] = Xv[ge][gl] + op_w + {7'b0, op_s_w};

                    wire        sgn_w = q_base[ge][gl][7];
                    wire [6:0]  lo_w  = sgn_w ? ~q_base[ge][gl][6:0] : q_base[ge][gl][6:0];
                    wire        sat_w = (lo_w[6:5] != 2'b00);
                    wire [4:0]  mp_w  = sat_w ? 5'd31 : lo_w[4:0];
                    // Forced to key 63 (the max possible value -- always
                    // loses ties, see merge_top2's comment) while this
                    // window's q_sgn is forced separately to 0.
                    assign q_key[ge][gl] = ((state == S_LOAD) && !overlap_cycle) ? 6'd63 : {mp_w, sgn_w};
                    assign q_sgn[ge][gl] = ((state == S_LOAD) && !overlap_cycle) ? 1'b0  : sgn_w;
                end else begin: LANE_NORMAL
                    assign x_base[ge][gl] = Xv[ge][gl] + op_w + {7'b0, op_s_w};

                    wire        sgn_w = q_base[ge][gl][7];
                    wire [6:0]  lo_w  = sgn_w ? ~q_base[ge][gl][6:0] : q_base[ge][gl][6:0];
                    wire        sat_w = (lo_w[6:5] != 2'b00);
                    wire [4:0]  mp_w  = sat_w ? 5'd31 : lo_w[4:0];
                    assign q_key[ge][gl] = {mp_w, sgn_w};
                    assign q_sgn[ge][gl] = sgn_w;
                end
            end
        end
    endgenerate

    wire [6:0] qneg_vec [0:15];
    genvar qgl;
    generate
        for (qgl = 0; qgl < 16; qgl = qgl + 1) begin: QNEG
            assign qneg_vec[qgl] = {q_sgn[6][qgl], q_sgn[5][qgl], q_sgn[4][qgl], q_sgn[3][qgl],
                                     q_sgn[2][qgl], q_sgn[1][qgl], q_sgn[0][qgl]};
        end
    endgenerate

    wire [19:0] c2v_new_lane   [0:15];
    wire [34:0] rmag_flat_lane [0:15];
    wire  [6:0] rneg_lane      [0:15];

    genvar glane;
    generate
        for (glane = 0; glane < 16; glane = glane + 1) begin: CNU_INST
            CNU_lane u_cnu (
                .key0(q_key[0][glane]), .key1(q_key[1][glane]),
                .key2(q_key[2][glane]), .key3(q_key[3][glane]),
                .key4(q_key[4][glane]), .key5(q_key[5][glane]),
                .key6(q_key[6][glane]),
                .qneg(qneg_vec[glane]),
                .c2v_new(c2v_new_lane[glane]),
                .rmag_flat(rmag_flat_lane[glane]),
                .rneg(rneg_lane[glane])
            );
        end
    endgenerate

    //============================================================
    // Fast path for lane 12 / the overlap cycle.
    // cap_desc captures c2v_new_lane[13] (the 6-real-edge, slot-7-
    // forced-loss result) one cycle before it is needed, so the
    // overlap cycle only has to combine it with in_data through a
    // short chain instead of a live 7-edge tree.
    //============================================================

    reg [19:0] cap_desc;
    always @(posedge clk) begin
        if (state == S_LOAD) cap_desc <= c2v_new_lane[13];
    end

    // S1's table, duplicated at top level for the fast path's norm(|in_data|).
    function [4:0] norm5_top;
        input [4:0] magnitude;
        begin
            case (magnitude)
                5'd0:  norm5_top = 5'd0;
                5'd1:  norm5_top = 5'd1;
                5'd2:  norm5_top = 5'd2;
                5'd3:  norm5_top = 5'd2;
                5'd4:  norm5_top = 5'd3;
                5'd5:  norm5_top = 5'd4;
                5'd6:  norm5_top = 5'd5;
                5'd7:  norm5_top = 5'd5;
                5'd8:  norm5_top = 5'd6;
                5'd9:  norm5_top = 5'd7;
                5'd10: norm5_top = 5'd8;
                5'd11: norm5_top = 5'd8;
                5'd12: norm5_top = 5'd9;
                5'd13: norm5_top = 5'd10;
                5'd14: norm5_top = 5'd11;
                5'd15: norm5_top = 5'd11;
                5'd16: norm5_top = 5'd12;
                5'd17: norm5_top = 5'd13;
                5'd18: norm5_top = 5'd14;
                5'd19: norm5_top = 5'd14;
                5'd20: norm5_top = 5'd15;
                5'd21: norm5_top = 5'd16;
                5'd22: norm5_top = 5'd17;
                5'd23: norm5_top = 5'd17;
                5'd24: norm5_top = 5'd18;
                5'd25: norm5_top = 5'd19;
                5'd26: norm5_top = 5'd20;
                5'd27: norm5_top = 5'd20;
                5'd28: norm5_top = 5'd21;
                5'd29: norm5_top = 5'd22;
                5'd30: norm5_top = 5'd23;
                default: norm5_top = 5'd23; // 31
            endcase
        end
    endfunction

    wire        xs    = in_data[5];
    wire [4:0]  x_mag = in_data[5] ? (~in_data[4:0] + 5'd1) : in_data[4:0];
    wire [4:0]  nx    = norm5_top(x_mag);

    wire [4:0] cn1   = cap_desc[19:15];
    wire [4:0] cn2   = cap_desc[14:10];
    wire [2:0] cidx  = cap_desc[9:7];
    wire [6:0] crneg = cap_desc[6:0];

    // v6-2 (B1): thr24(c) = min{m : norm(m) >= c}. Since norm is monotonic
    // non-decreasing, "nx < c" is equivalent to "|in_data| < thr24(c)".
    // thr1/thr2 depend only on cap_desc (an FF, stable well before in_data
    // arrives), so in_data is compared directly against +-thr instead of
    // waiting for its own abs+norm (nx) to resolve first.
    function [4:0] thr24;
        input [4:0] c;
        begin
            case (c)
                5'd0:  thr24 = 5'd0;
                5'd1:  thr24 = 5'd1;
                5'd2:  thr24 = 5'd2;
                5'd3:  thr24 = 5'd4;
                5'd4:  thr24 = 5'd5;
                5'd5:  thr24 = 5'd6;
                5'd6:  thr24 = 5'd8;
                5'd7:  thr24 = 5'd9;
                5'd8:  thr24 = 5'd10;
                5'd9:  thr24 = 5'd12;
                5'd10: thr24 = 5'd13;
                5'd11: thr24 = 5'd14;
                5'd12: thr24 = 5'd16;
                5'd13: thr24 = 5'd17;
                5'd14: thr24 = 5'd18;
                5'd15: thr24 = 5'd20;
                5'd16: thr24 = 5'd21;
                5'd17: thr24 = 5'd22;
                5'd18: thr24 = 5'd24;
                5'd19: thr24 = 5'd25;
                5'd20: thr24 = 5'd26;
                5'd21: thr24 = 5'd28;
                5'd22: thr24 = 5'd29;
                default: thr24 = 5'd30; // c == 23
            endcase
        end
    endfunction

    wire [4:0]        thr1     = thr24(cn1);
    wire [4:0]        thr2     = thr24(cn2);
    wire signed [6:0] in_data_s = {in_data[5], in_data};
    wire signed [6:0] neg_thr1  = -{2'b00, thr1};
    wire signed [6:0] neg_thr2  = -{2'b00, thr2};
    wire lt1 = in_data[5] ? (in_data_s > neg_thr1) : (in_data[4:0] < thr1);
    wire lt2 = in_data[5] ? (in_data_s > neg_thr2) : (in_data[4:0] < thr2);

    // v5 (A): the two comparisons against in_data's magnitude are shared.
    // min(select(cn1, cn2), nx) == select(min(cn1, nx), min(cn2, nx)), so
    // each edge only picks between p1/p2 instead of owning a comparator.
    wire [4:0] p1 = lt1 ? nx : cn1;
    wire [4:0] p2 = lt2 ? nx : cn2;

    // Per-edge fast-path r_new for lane 12 (used by newv below): edge 6's
    // own value never depends on itself. v6-2's cidx!=6 simplification:
    // its capture-forced key (63, the maximum possible value) can never
    // win a tie (see CNU_lane's merge_top2 comment), so cidx is never 6 --
    // rmag_fast[6] is simply cn1 (crneg[6] unchanged), no mux needed. The
    // other 6 edges must now also consider the fresh in_data.
    wire [4:0] rmag_fast [0:6];
    wire       rneg_fast [0:6];
    genvar fe;
    generate
        for (fe = 0; fe < 7; fe = fe + 1) begin: FASTEDGE
            if (fe == 6) begin: FAST_SELF
                assign rmag_fast[fe] = cn1;
                assign rneg_fast[fe] = crneg[fe];
            end else begin: FAST_OTHER
                assign rmag_fast[fe] = (cidx == fe) ? p2 : p1;
                assign rneg_fast[fe] = crneg[fe] ^ xs;
            end
        end
    endgenerate

    // New FIFO descriptor for lane 12: the true top-2/idx among all 7
    // edges now that in_data (edge 6) is a real candidate.
    wire [4:0] fast_m1   = p1;
    wire [4:0] fast_m2   = lt1 ? cn1 : p2;
    wire [2:0] fast_idx  = lt1 ? 3'd6 : cidx;
    wire [6:0] fast_rneg = {rneg_fast[6], rneg_fast[5], rneg_fast[4], rneg_fast[3],
                            rneg_fast[2], rneg_fast[1], rneg_fast[0]};
    wire [19:0] fast_desc = {fast_m1, fast_m2, fast_idx, fast_rneg};

    wire signed [7:0] newv [0:6][0:15];

    genvar ne, nl;
    generate
        for (ne = 0; ne < 7; ne = ne + 1) begin: NEWV_EDGE
            for (nl = 0; nl < 16; nl = nl + 1) begin: NEWV_LANE
                if (nl == 12) begin: LANE12_NEWV
                    wire [4:0] rmag_e   = rmag_flat_lane[nl][ne*5 +: 5];
                    wire       rneg_e   = rneg_lane[nl][ne];
                    wire [4:0] rmag_sel = overlap_cycle ? rmag_fast[ne] : rmag_e;
                    wire       rneg_sel = overlap_cycle ? rneg_fast[ne] : rneg_e;
                    wire [7:0] rop_e    = {3'b0, rmag_sel} ^ {8{rneg_sel}};
                    assign newv[ne][nl] = x_base[ne][nl] + rop_e + {7'b0, rneg_sel};
                end else begin: NORMAL_NEWV
                    wire [4:0] rmag_e = rmag_flat_lane[nl][ne*5 +: 5];
                    wire       rneg_e = rneg_lane[nl][ne];
                    wire [7:0] rop_e  = {3'b0, rmag_e} ^ {8{rneg_e}};
                    assign newv[ne][nl] = x_base[ne][nl] + rop_e + {7'b0, rneg_e};
                end
            end
        end
    endgenerate

    //============================================================
    // c2v FIFO: 4 stages, unconditional shift every cycle. While idle,
    // the tail clears in full -- idx and rneg are read directly by
    // ro_neg's bit-select, not through a magnitude-based mux, so leaving
    // them uninitialized would poison q_base/x_base through op_s_w (see
    // arch_notes.md 4.2's S3 note: this is a 4-state-simulation-only
    // hazard, but it is real for RTL/gate-level sim).
    //============================================================

    always @(posedge clk) begin
        for (rr = 0; rr < 16; rr = rr + 1) begin
            c2v_fifo[0][rr] <= c2v_fifo[1][rr];
            c2v_fifo[1][rr] <= c2v_fifo[2][rr];
            c2v_fifo[2][rr] <= c2v_fifo[3][rr];
            if (run_now) begin
                c2v_fifo[3][rr] <= (rr == 12 && overlap_cycle) ? fast_desc : c2v_new_lane[rr];
            end else begin
                c2v_fifo[3][rr] <= 20'd0;
            end
        end
    end

    //============================================================
    // Bank A / B storage per column.
    //
    // S8 (A): every column does something every cycle -- no dedicated
    // "hold" leg anywhere. Priority per column: (1) commit -> the
    // touch-based schedule below (self-rotate / new-write / an explicit
    // rot1 for an absent layer that needs compensation); (2) this
    // column's own LOAD injection window (rotate-by-1, overwriting one
    // position with in_data); (3) otherwise, plain rotate-by-1 -- this
    // covers IDLE, LOAD cycles that are not this column's turn, all of
    // OUTPUT, and the cycle where RUN decides to stop.
    // Injection points inj[c] = (-tf[c]-2) mod 16 = [9,0,4,12,1,2,5,11]
    // (one earlier than v3's, to compensate for every column rotating
    // through the other columns' LOAD windows too).
    //
    // S6 (B): only columns 1 and 2 keep an explicit hold (they are the
    // only two that hit an absent layer while their B value is still
    // between first and last touch). Every other column's B write is
    // an unconditional select among its real sources, never gated.
    //============================================================

    wire signed [7:0] inj_val = {{2{in_data[5]}}, in_data};

    // ---- column 0 : inject 9, read 11, edge = layer-1
    always @(posedge clk) begin
        if (commit) begin
            if (nl_r1) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[0][rr] <= A[0][(rr + 8) % 16];
            end else if (nl_r3) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[0][rr] <= newv[2][(rr + 8) % 16];
            end
            // layer0 absent (y=0, reuses the existing delta-0/hold leg);
            // layer2 self-rotate delta 0 (same leg).
        end else if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd0) && (io_cnt != 7'd127)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[0][rr] <= (rr == 9) ? inj_val : A[0][(rr + 1) % 16];
        end else begin
            for (rr = 0; rr < 16; rr = rr + 1) A[0][rr] <= A[0][(rr + 1) % 16];
        end
    end

    always @(posedge clk) begin
        for (rr = 0; rr < 16; rr = rr + 1)
            B[0][rr] <= nl_r1 ? newv[0][(rr + 8) % 16] : newv[1][rr];
    end

    // ---- column 1 : inject 0, read 2, edge = 0
    always @(posedge clk) begin
        if (commit) begin
            if (nl_r0 || nl_r2) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[1][rr] <= A[1][(rr + 4) % 16];
            end else if (nl_r3) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[1][rr] <= newv[0][(rr + 8) % 16];
            end
            // layer1 absent (y=0, reuses the existing delta-0/hold leg)
        end else if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd1) && (io_cnt != 7'd127)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[1][rr] <= (rr == 0) ? inj_val : A[1][(rr + 1) % 16];
        end else begin
            for (rr = 0; rr < 16; rr = rr + 1) A[1][rr] <= A[1][(rr + 1) % 16];
        end
    end

    always @(posedge clk) begin
        if (run_now) begin
            if (nl_r0 || nl_r2) begin
                for (rr = 0; rr < 16; rr = rr + 1) B[1][rr] <= newv[0][(rr + 4) % 16];
            end
        end
    end

    // ---- column 2 : inject 4, read 6, edge = 1
    always @(posedge clk) begin
        if (commit) begin
            if (nl_r0) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[2][rr] <= A[2][(rr + 4) % 16];
            end else if (nl_r1) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[2][rr] <= A[2][(rr + 12) % 16];
            end else if (nl_r2) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[2][rr] <= A[2][(rr + 1) % 16]; // absent, y=1
            end else if (nl_r3) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[2][rr] <= newv[1][(rr + 15) % 16];
            end
        end else if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd2) && (io_cnt != 7'd127)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[2][rr] <= (rr == 4) ? inj_val : A[2][(rr + 1) % 16];
        end else begin
            for (rr = 0; rr < 16; rr = rr + 1) A[2][rr] <= A[2][(rr + 1) % 16];
        end
    end

    always @(posedge clk) begin
        if (run_now) begin
            if (nl_r0) begin
                for (rr = 0; rr < 16; rr = rr + 1) B[2][rr] <= newv[1][(rr + 4) % 16];
            end else if (nl_r1) begin
                for (rr = 0; rr < 16; rr = rr + 1) B[2][rr] <= newv[1][(rr + 13) % 16];
            end
        end
    end

    // ---- column 3 : inject 12, read 14, edge = 2
    always @(posedge clk) begin
        if (commit) begin
            if (nl_r0) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[3][rr] <= A[3][(rr + 8) % 16];
            end else if (nl_r1) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[3][rr] <= A[3][(rr + 1) % 16];
            end else if (nl_r2) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[3][rr] <= newv[2][(rr + 6) % 16];
            end else if (nl_r3) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[3][rr] <= A[3][(rr + 1) % 16]; // absent, y=1
            end
        end else if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd3) && (io_cnt != 7'd127)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[3][rr] <= (rr == 12) ? inj_val : A[3][(rr + 1) % 16];
        end else begin
            for (rr = 0; rr < 16; rr = rr + 1) A[3][rr] <= A[3][(rr + 1) % 16];
        end
    end

    always @(posedge clk) begin
        for (rr = 0; rr < 16; rr = rr + 1)
            B[3][rr] <= nl_r0 ? newv[2][(rr + 8) % 16] : newv[2][(rr + 1) % 16];
    end

    // ---- column 4 : inject 1, read 3, edge = 3, present every layer
    always @(posedge clk) begin
        if (commit) begin
            if (nl_r0 || nl_r1) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[4][rr] <= A[4][(rr + 5) % 16];
            end else if (nl_r2) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[4][rr] <= A[4][(rr + 13) % 16];
            end else if (nl_r3) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[4][rr] <= newv[3][(rr + 9) % 16];
            end
        end else if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd4) && (io_cnt != 7'd127)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[4][rr] <= (rr == 1) ? inj_val : A[4][(rr + 1) % 16];
        end else begin
            for (rr = 0; rr < 16; rr = rr + 1) A[4][rr] <= A[4][(rr + 1) % 16];
        end
    end

    always @(posedge clk) begin
        for (rr = 0; rr < 16; rr = rr + 1)
            B[4][rr] <= nl_r2 ? newv[3][(rr + 13) % 16] : newv[3][(rr + 5) % 16];
    end

    // ---- column 5 : inject 2, read 4, edge = 4, present every layer
    always @(posedge clk) begin
        if (commit) begin
            if (nl_r0 || nl_r2) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[5][rr] <= A[5][(rr + 1) % 16];
            end else if (nl_r1) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[5][rr] <= A[5][(rr + 2) % 16];
            end else if (nl_r3) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[5][rr] <= newv[4][(rr + 12) % 16];
            end
        end else if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd5) && (io_cnt != 7'd127)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[5][rr] <= (rr == 2) ? inj_val : A[5][(rr + 1) % 16];
        end else begin
            for (rr = 0; rr < 16; rr = rr + 1) A[5][rr] <= A[5][(rr + 1) % 16];
        end
    end

    always @(posedge clk) begin
        for (rr = 0; rr < 16; rr = rr + 1)
            B[5][rr] <= nl_r1 ? newv[4][(rr + 2) % 16] : newv[4][(rr + 1) % 16];
    end

    // ---- column 6 : inject 5, read 7, edge = 5, present every layer
    always @(posedge clk) begin
        if (commit) begin
            if (nl_r0) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[6][rr] <= A[6][(rr + 3) % 16];
            end else if (nl_r1 || nl_r2) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[6][rr] <= A[6][(rr + 14) % 16];
            end else if (nl_r3) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[6][rr] <= newv[5][(rr + 1) % 16];
            end
        end else if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd6) && (io_cnt != 7'd127)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[6][rr] <= (rr == 5) ? inj_val : A[6][(rr + 1) % 16];
        end else begin
            for (rr = 0; rr < 16; rr = rr + 1) A[6][rr] <= A[6][(rr + 1) % 16];
        end
    end

    always @(posedge clk) begin
        for (rr = 0; rr < 16; rr = rr + 1)
            B[6][rr] <= nl_r0 ? newv[5][(rr + 3) % 16] : newv[5][(rr + 14) % 16];
    end

    // ---- column 7 : inject 11, read 13, edge = 6, present every layer.
    // Only j=0..14 (io_cnt 112..126) rotate+inject; j=15 (Lch127) is
    // folded into the overlap cycle instead. A_7[6]'s self-rotate source
    // (normally A_7[12]) is replaced by in_data directly during that cycle.
    always @(posedge clk) begin
        if (commit) begin
            if (nl_r0) begin
                for (rr = 0; rr < 16; rr = rr + 1)
                    A[7][rr] <= (rr == 6 && overlap_cycle) ? inj_val : A[7][(rr + 6) % 16];
            end else if (nl_r2) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[7][rr] <= A[7][(rr + 10) % 16];
            end else if (nl_r3) begin
                for (rr = 0; rr < 16; rr = rr + 1) A[7][rr] <= newv[6][rr];
            end
            // layer1 self-rotate delta 0 (reuses the existing hold leg)
        end else if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd7) && (io_cnt != 7'd127)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[7][rr] <= (rr == 11) ? inj_val : A[7][(rr + 1) % 16];
        end else begin
            for (rr = 0; rr < 16; rr = rr + 1) A[7][rr] <= A[7][(rr + 1) % 16];
        end
    end

    always @(posedge clk) begin
        for (rr = 0; rr < 16; rr = rr + 1)
            B[7][rr] <= nl_r0 ? newv[6][(rr + 6) % 16] : (nl_r1 ? newv[6][rr] : newv[6][(rr + 10) % 16]);
    end

    //============================================================
    // Syndrome network: fixed XOR wiring on A's sign bits.
    // rot_amt[n][c] = (BG[n][c] - tf[c]) mod 16, tf = [5,14,10,2,13,12,9,3].
    // Unaffected by S8: A's alignment at every iteration boundary is the
    // same as v1-v3 (t_f), only the between-boundary rotation schedule
    // changed.
    //============================================================

    function [15:0] rotate_read16;
        input [15:0] value;
        input [3:0]  shift;
        reg [31:0] doubled;
        begin
            doubled = {value, value};
            rotate_read16 = doubled >> shift;
        end
    endfunction

    wire [15:0] hard [0:7];
    genvar hc, hr;
    generate
        for (hc = 0; hc < 8; hc = hc + 1) begin: HARD_COL
            for (hr = 0; hr < 16; hr = hr + 1) begin: HARD_ROW
                assign hard[hc][hr] = A[hc][hr][7];
            end
        end
    endgenerate

    wire [15:0] syn0, syn1, syn2, syn3;

    assign syn0 = hard[1] ^ hard[2] ^ hard[3] ^ hard[4] ^ hard[5] ^ hard[6] ^ hard[7];

    assign syn1 = hard[0]
                ^ rotate_read16(hard[2], 4'd4)
                ^ rotate_read16(hard[3], 4'd8)
                ^ rotate_read16(hard[4], 4'd5)
                ^ rotate_read16(hard[5], 4'd1)
                ^ rotate_read16(hard[6], 4'd3)
                ^ rotate_read16(hard[7], 4'd6);

    assign syn2 = rotate_read16(hard[0], 4'd11)
                ^ rotate_read16(hard[1], 4'd7)
                ^ rotate_read16(hard[3], 4'd12)
                ^ rotate_read16(hard[4], 4'd13)
                ^ rotate_read16(hard[5], 4'd6)
                ^ rotate_read16(hard[6], 4'd4)
                ^ rotate_read16(hard[7], 4'd9);

    assign syn3 = rotate_read16(hard[0], 4'd2)
                ^ rotate_read16(hard[1], 4'd2)
                ^ rotate_read16(hard[2], 4'd11)
                ^ rotate_read16(hard[4], 4'd1)
                ^ rotate_read16(hard[5], 4'd14)
                ^ rotate_read16(hard[6], 4'd9)
                ^ rotate_read16(hard[7], 4'd10);

    assign syndrome_nonzero = |(syn0 | syn1 | syn2 | syn3);

    //============================================================
    // Output: rotate-by-1 read, fixed point per column, 8:1 column mux.
    // pr = [11,2,6,14,3,4,7,13] for columns 0..7 -- unchanged by S8
    // (each column has already rotated 16c times before its own output
    // window starts, i.e. a full cycle, so the fixed read point is the
    // same as v1-v3).
    //============================================================

    wire [2:0] eff_col = stop_now ? 3'd0 : io_cnt[6:4];

    reg signed [7:0] out_data_comb;
    always @(*) begin
        case (eff_col)
            3'd0: out_data_comb = A[0][11];
            3'd1: out_data_comb = A[1][2];
            3'd2: out_data_comb = A[2][6];
            3'd3: out_data_comb = A[3][14];
            3'd4: out_data_comb = A[4][3];
            3'd5: out_data_comb = A[5][4];
            3'd6: out_data_comb = A[6][7];
            default: out_data_comb = A[7][13];
        endcase
    end

    reg warn_reg;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            out_valid <= 1'b0;
            out_data  <= 8'd0;
            out_warn  <= 1'b0;
            warn_reg  <= 1'b0;
        end else if (stop_now) begin
            warn_reg  <= at_limit && syndrome_nonzero;
            out_valid <= 1'b1;
            out_warn  <= at_limit && syndrome_nonzero;
            out_data  <= out_data_comb;
        end else if (state == S_OUT) begin
            out_valid <= 1'b1;
            out_warn  <= warn_reg;
            out_data  <= out_data_comb;
        end else begin
            out_valid <= 1'b0;
            out_data  <= 8'd0;
            out_warn  <= 1'b0;
        end
    end

endmodule
