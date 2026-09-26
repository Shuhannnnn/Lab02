// LDPC.v -- Lab02 v1 architecture.
// 16 CNU lanes (one full layer per cycle), rotated storage in bank A/B,
// touch-based (first/mid/last) bank control shared by both modes,
// 4-stage c2v FIFO, rotate-by-1 input load / output.
// See arch_notes.md section 3 for the derivation of every constant below.

module CNU_lane (
    input  signed [5:0] q0,
    input  signed [5:0] q1,
    input  signed [5:0] q2,
    input  signed [5:0] q3,
    input  signed [5:0] q4,
    input  signed [5:0] q5,
    input  signed [5:0] q6,

    output reg [19:0] c2v_new,
    output reg [41:0] r_new_flat
);

    reg [12:0] leaf [0:6];
    reg [12:0] pair01;
    reg [12:0] pair23;
    reg [12:0] pair45;
    reg [12:0] group03;
    reg [12:0] group46;
    reg [12:0] root_pair;

    reg [4:0] mag0;
    reg [4:0] mag1;
    reg [4:0] mag2;
    reg [4:0] mag3;
    reg [4:0] mag4;
    reg [4:0] mag5;
    reg [4:0] mag6;
    reg [6:0] q_neg;
    reg [6:0] r_neg;
    reg       sign_parity;
    reg [4:0] norm_min1;
    reg [4:0] norm_min2;

    function [4:0] abs6;
        input signed [5:0] value;
        begin
            if (value[5])
                abs6 = (~value[4:0]) + 5'd1;
            else
                abs6 = value[4:0];
        end
    endfunction

    function [12:0] merge_top2;
        input [12:0] left_value;
        input [12:0] right_value;
        reg [4:0] next_min2;
        begin
            // The tree always places lower edge indices in the left subtree.
            // Selecting left on a tie therefore implements the required tie rule.
            if (left_value[12:8] <= right_value[12:8]) begin
                if (left_value[4:0] < right_value[12:8])
                    next_min2 = left_value[4:0];
                else
                    next_min2 = right_value[12:8];
                merge_top2 = {left_value[12:5], next_min2};
            end else begin
                if (right_value[4:0] < left_value[12:8])
                    next_min2 = right_value[4:0];
                else
                    next_min2 = left_value[12:8];
                merge_top2 = {right_value[12:5], next_min2};
            end
        end
    endfunction

    function [4:0] norm5;
        input [4:0] magnitude;
        reg [6:0] wide_value;
        begin
            wide_value = ({2'b00, magnitude} << 1)
                       +  {2'b00, magnitude} + 7'd2;
            norm5 = wide_value[6:2];
        end
    endfunction

    function [5:0] make_r;
        input [4:0] magnitude;
        input       negative;
        begin
            if (magnitude == 5'd0)
                make_r = 6'd0;
            else if (negative)
                make_r = (~{1'b0, magnitude}) + 6'd1;
            else
                make_r = {1'b0, magnitude};
        end
    endfunction

    always @(*) begin
        mag0 = abs6(q0);
        mag1 = abs6(q1);
        mag2 = abs6(q2);
        mag3 = abs6(q3);
        mag4 = abs6(q4);
        mag5 = abs6(q5);
        mag6 = abs6(q6);

        leaf[0] = {mag0, 3'd0, 5'd31};
        leaf[1] = {mag1, 3'd1, 5'd31};
        leaf[2] = {mag2, 3'd2, 5'd31};
        leaf[3] = {mag3, 3'd3, 5'd31};
        leaf[4] = {mag4, 3'd4, 5'd31};
        leaf[5] = {mag5, 3'd5, 5'd31};
        leaf[6] = {mag6, 3'd6, 5'd31};

        pair01   = merge_top2(leaf[0], leaf[1]);
        pair23   = merge_top2(leaf[2], leaf[3]);
        pair45   = merge_top2(leaf[4], leaf[5]);
        group03  = merge_top2(pair01, pair23);
        group46  = merge_top2(pair45, leaf[6]);
        root_pair = merge_top2(group03, group46);

        norm_min1 = norm5(root_pair[12:8]);
        norm_min2 = norm5(root_pair[4:0]);

        q_neg = {q6[5], q5[5], q4[5], q3[5], q2[5], q1[5], q0[5]};
        sign_parity = ^q_neg;
        r_neg = q_neg ^ {7{sign_parity}};

        c2v_new = {norm_min1, norm_min2, root_pair[7:5], r_neg};

        r_new_flat[5:0] = make_r(
            (root_pair[7:5] == 3'd0) ? norm_min2 : norm_min1,
            r_neg[0]);
        r_new_flat[11:6] = make_r(
            (root_pair[7:5] == 3'd1) ? norm_min2 : norm_min1,
            r_neg[1]);
        r_new_flat[17:12] = make_r(
            (root_pair[7:5] == 3'd2) ? norm_min2 : norm_min1,
            r_neg[2]);
        r_new_flat[23:18] = make_r(
            (root_pair[7:5] == 3'd3) ? norm_min2 : norm_min1,
            r_neg[3]);
        r_new_flat[29:24] = make_r(
            (root_pair[7:5] == 3'd4) ? norm_min2 : norm_min1,
            r_neg[4]);
        r_new_flat[35:30] = make_r(
            (root_pair[7:5] == 3'd5) ? norm_min2 : norm_min1,
            r_neg[5]);
        r_new_flat[41:36] = make_r(
            (root_pair[7:5] == 3'd6) ? norm_min2 : norm_min1,
            r_neg[6]);
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
    wire commit   = (state == S_RUN) && !stop_now;

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
        end else if ((state == S_RUN) && !stop_now) begin
            layer <= (layer == 2'd3) ? 2'd0 : layer + 2'd1;
            if (layer == 2'd3) iter <= iter + 4'd1;
        end
    end

    //============================================================
    // Per-layer edge table: edge e (0..6) is normally column e+1;
    // column 0 borrows edge (layer-1) at layer 1/2/3.
    //============================================================

    reg [2:0] col_sel    [0:6];
    reg       touch_first [0:6];
    integer   ti;

    always @(*) begin
        case (layer)
            2'd0: begin
                col_sel[0] = 3'd1; col_sel[1] = 3'd2; col_sel[2] = 3'd3;
                col_sel[3] = 3'd4; col_sel[4] = 3'd5; col_sel[5] = 3'd6;
                col_sel[6] = 3'd7;
                for (ti = 0; ti < 7; ti = ti + 1) touch_first[ti] = 1'b1;
            end
            2'd1: begin
                col_sel[0] = 3'd0; col_sel[1] = 3'd2; col_sel[2] = 3'd3;
                col_sel[3] = 3'd4; col_sel[4] = 3'd5; col_sel[5] = 3'd6;
                col_sel[6] = 3'd7;
                touch_first[0] = 1'b1;
                for (ti = 1; ti < 7; ti = ti + 1) touch_first[ti] = 1'b0;
            end
            2'd2: begin
                col_sel[0] = 3'd1; col_sel[1] = 3'd0; col_sel[2] = 3'd3;
                col_sel[3] = 3'd4; col_sel[4] = 3'd5; col_sel[5] = 3'd6;
                col_sel[6] = 3'd7;
                for (ti = 0; ti < 7; ti = ti + 1) touch_first[ti] = 1'b0;
            end
            default: begin // layer 3
                col_sel[0] = 3'd1; col_sel[1] = 3'd2; col_sel[2] = 3'd0;
                col_sel[3] = 3'd4; col_sel[4] = 3'd5; col_sel[5] = 3'd6;
                col_sel[6] = 3'd7;
                for (ti = 0; ti < 7; ti = ti + 1) touch_first[ti] = 1'b0;
            end
        endcase
    end

    //============================================================
    // Edge datapath (7 edges x 16 lanes)
    //============================================================

    wire signed [7:0] Xv     [0:6][0:15];
    wire signed [7:0] Qv     [0:6][0:15];
    wire        [4:0] ro_mag [0:6][0:15];
    wire               ro_neg [0:6][0:15];
    wire signed [7:0] ro     [0:6][0:15];
    wire signed [7:0] q_base [0:6][0:15];
    wire signed [7:0] x_base [0:6][0:15];
    wire signed [5:0] q_msg  [0:6][0:15];

    function [5:0] clip6;
        input signed [7:0] value;
        begin
            if (value > 8'sd31)
                clip6 = 6'sd31;
            else if (value < -8'sd31)
                clip6 = -6'sd31;
            else
                clip6 = value[5:0];
        end
    endfunction

    genvar ge, gl;
    generate
        for (ge = 0; ge < 7; ge = ge + 1) begin: EDGE
            for (gl = 0; gl < 16; gl = gl + 1) begin: LANE
                assign Xv[ge][gl] = touch_first[ge] ? A[col_sel[ge]][gl]
                                                     : B[col_sel[ge]][gl];
                assign Qv[ge][gl] = (!mode_reg || touch_first[ge]) ? A[col_sel[ge]][gl]
                                                                    : B[col_sel[ge]][gl];
                assign ro_mag[ge][gl] = (c2v_fifo[0][gl][9:7] == ge[2:0])
                                        ? c2v_fifo[0][gl][14:10]
                                        : c2v_fifo[0][gl][19:15];
                assign ro_neg[ge][gl] = c2v_fifo[0][gl][ge];
                assign ro[ge][gl] = ro_neg[ge][gl] ? -$signed({3'b0, ro_mag[ge][gl]})
                                                    :  $signed({3'b0, ro_mag[ge][gl]});
                assign q_base[ge][gl] = Qv[ge][gl] - ro[ge][gl];
                assign x_base[ge][gl] = Xv[ge][gl] - ro[ge][gl];
                assign q_msg[ge][gl]  = clip6(q_base[ge][gl]);
            end
        end
    endgenerate

    wire [19:0] c2v_new_lane    [0:15];
    wire [41:0] r_new_flat_lane [0:15];

    genvar glane;
    generate
        for (glane = 0; glane < 16; glane = glane + 1) begin: CNU_INST
            CNU_lane u_cnu (
                .q0(q_msg[0][glane]), .q1(q_msg[1][glane]),
                .q2(q_msg[2][glane]), .q3(q_msg[3][glane]),
                .q4(q_msg[4][glane]), .q5(q_msg[5][glane]),
                .q6(q_msg[6][glane]),
                .c2v_new(c2v_new_lane[glane]),
                .r_new_flat(r_new_flat_lane[glane])
            );
        end
    endgenerate

    wire signed [7:0] newv [0:6][0:15];

    genvar ne, nl;
    generate
        for (ne = 0; ne < 7; ne = ne + 1) begin: NEWV_EDGE
            for (nl = 0; nl < 16; nl = nl + 1) begin: NEWV_LANE
                wire signed [5:0] r_new_e = r_new_flat_lane[nl][ne*6 +: 6];
                assign newv[ne][nl] = x_base[ne][nl] + {{2{r_new_e[5]}}, r_new_e};
            end
        end
    endgenerate

    //============================================================
    // c2v FIFO: 4 stages, unconditional shift every cycle.
    //============================================================

    always @(posedge clk) begin
        for (rr = 0; rr < 16; rr = rr + 1) begin
            c2v_fifo[0][rr] <= c2v_fifo[1][rr];
            c2v_fifo[1][rr] <= c2v_fifo[2][rr];
            c2v_fifo[2][rr] <= c2v_fifo[3][rr];
            c2v_fifo[3][rr] <= (state == S_RUN) ? c2v_new_lane[rr] : 20'd0;
        end
    end

    //============================================================
    // Bank A / B storage per column.
    // Load/output use rotate-by-1 with one fixed injection/read point.
    // RUN uses touch-based self-rotate (A) / new-write (A on last touch,
    // B on first/mid touch), driven by the per-layer edge table above.
    //============================================================

    wire signed [7:0] inj_val = {{2{in_data[5]}}, in_data};

    // ---- column 0 : inject 10, read 11, edge = layer-1
    always @(posedge clk) begin
        if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd0)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[0][rr] <= (rr == 10) ? inj_val : A[0][(rr + 1) % 16];
        end else if (((state == S_OUT) && (io_cnt[6:4] == 3'd0)) || stop_now) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[0][rr] <= A[0][(rr + 1) % 16];
        end else if (commit) begin
            case (layer)
                2'd1: for (rr = 0; rr < 16; rr = rr + 1) A[0][rr] <= A[0][(rr + 8) % 16];
                2'd3: for (rr = 0; rr < 16; rr = rr + 1) A[0][rr] <= newv[2][(rr + 8) % 16];
                default: ; // layer0 absent; layer2 self-rotate delta 0 (hold)
            endcase
        end
    end

    always @(posedge clk) begin
        if (state == S_RUN) begin
            case (layer)
                2'd1: for (rr = 0; rr < 16; rr = rr + 1) B[0][rr] <= newv[0][(rr + 8) % 16];
                2'd2: for (rr = 0; rr < 16; rr = rr + 1) B[0][rr] <= newv[1][rr];
                default: ;
            endcase
        end
    end

    // ---- column 1 : inject 1, read 2, edge = 0
    always @(posedge clk) begin
        if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd1)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[1][rr] <= (rr == 1) ? inj_val : A[1][(rr + 1) % 16];
        end else if ((state == S_OUT) && (io_cnt[6:4] == 3'd1)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[1][rr] <= A[1][(rr + 1) % 16];
        end else if (commit) begin
            case (layer)
                2'd0: for (rr = 0; rr < 16; rr = rr + 1) A[1][rr] <= A[1][(rr + 4) % 16];
                2'd2: for (rr = 0; rr < 16; rr = rr + 1) A[1][rr] <= A[1][(rr + 4) % 16];
                2'd3: for (rr = 0; rr < 16; rr = rr + 1) A[1][rr] <= newv[0][(rr + 8) % 16];
                default: ; // layer1 absent
            endcase
        end
    end

    always @(posedge clk) begin
        if (state == S_RUN) begin
            case (layer)
                2'd0: for (rr = 0; rr < 16; rr = rr + 1) B[1][rr] <= newv[0][(rr + 4) % 16];
                2'd2: for (rr = 0; rr < 16; rr = rr + 1) B[1][rr] <= newv[0][(rr + 4) % 16];
                default: ;
            endcase
        end
    end

    // ---- column 2 : inject 5, read 6, edge = 1
    always @(posedge clk) begin
        if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd2)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[2][rr] <= (rr == 5) ? inj_val : A[2][(rr + 1) % 16];
        end else if ((state == S_OUT) && (io_cnt[6:4] == 3'd2)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[2][rr] <= A[2][(rr + 1) % 16];
        end else if (commit) begin
            case (layer)
                2'd0: for (rr = 0; rr < 16; rr = rr + 1) A[2][rr] <= A[2][(rr + 4) % 16];
                2'd1: for (rr = 0; rr < 16; rr = rr + 1) A[2][rr] <= A[2][(rr + 13) % 16];
                2'd3: for (rr = 0; rr < 16; rr = rr + 1) A[2][rr] <= newv[1][(rr + 15) % 16];
                default: ; // layer2 absent
            endcase
        end
    end

    always @(posedge clk) begin
        if (state == S_RUN) begin
            case (layer)
                2'd0: for (rr = 0; rr < 16; rr = rr + 1) B[2][rr] <= newv[1][(rr + 4) % 16];
                2'd1: for (rr = 0; rr < 16; rr = rr + 1) B[2][rr] <= newv[1][(rr + 13) % 16];
                default: ;
            endcase
        end
    end

    // ---- column 3 : inject 13, read 14, edge = 2
    always @(posedge clk) begin
        if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd3)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[3][rr] <= (rr == 13) ? inj_val : A[3][(rr + 1) % 16];
        end else if ((state == S_OUT) && (io_cnt[6:4] == 3'd3)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[3][rr] <= A[3][(rr + 1) % 16];
        end else if (commit) begin
            case (layer)
                2'd0: for (rr = 0; rr < 16; rr = rr + 1) A[3][rr] <= A[3][(rr + 8) % 16];
                2'd1: for (rr = 0; rr < 16; rr = rr + 1) A[3][rr] <= A[3][(rr + 1) % 16];
                2'd2: for (rr = 0; rr < 16; rr = rr + 1) A[3][rr] <= newv[2][(rr + 7) % 16];
                default: ; // layer3 absent
            endcase
        end
    end

    always @(posedge clk) begin
        if (state == S_RUN) begin
            case (layer)
                2'd0: for (rr = 0; rr < 16; rr = rr + 1) B[3][rr] <= newv[2][(rr + 8) % 16];
                2'd1: for (rr = 0; rr < 16; rr = rr + 1) B[3][rr] <= newv[2][(rr + 1) % 16];
                default: ;
            endcase
        end
    end

    // ---- column 4 : inject 2, read 3, edge = 3, present every layer
    always @(posedge clk) begin
        if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd4)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[4][rr] <= (rr == 2) ? inj_val : A[4][(rr + 1) % 16];
        end else if ((state == S_OUT) && (io_cnt[6:4] == 3'd4)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[4][rr] <= A[4][(rr + 1) % 16];
        end else if (commit) begin
            case (layer)
                2'd0: for (rr = 0; rr < 16; rr = rr + 1) A[4][rr] <= A[4][(rr + 5) % 16];
                2'd1: for (rr = 0; rr < 16; rr = rr + 1) A[4][rr] <= A[4][(rr + 5) % 16];
                2'd2: for (rr = 0; rr < 16; rr = rr + 1) A[4][rr] <= A[4][(rr + 13) % 16];
                2'd3: for (rr = 0; rr < 16; rr = rr + 1) A[4][rr] <= newv[3][(rr + 9) % 16];
            endcase
        end
    end

    always @(posedge clk) begin
        if (state == S_RUN) begin
            case (layer)
                2'd0: for (rr = 0; rr < 16; rr = rr + 1) B[4][rr] <= newv[3][(rr + 5) % 16];
                2'd1: for (rr = 0; rr < 16; rr = rr + 1) B[4][rr] <= newv[3][(rr + 5) % 16];
                2'd2: for (rr = 0; rr < 16; rr = rr + 1) B[4][rr] <= newv[3][(rr + 13) % 16];
                default: ;
            endcase
        end
    end

    // ---- column 5 : inject 3, read 4, edge = 4, present every layer
    always @(posedge clk) begin
        if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd5)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[5][rr] <= (rr == 3) ? inj_val : A[5][(rr + 1) % 16];
        end else if ((state == S_OUT) && (io_cnt[6:4] == 3'd5)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[5][rr] <= A[5][(rr + 1) % 16];
        end else if (commit) begin
            case (layer)
                2'd0: for (rr = 0; rr < 16; rr = rr + 1) A[5][rr] <= A[5][(rr + 1) % 16];
                2'd1: for (rr = 0; rr < 16; rr = rr + 1) A[5][rr] <= A[5][(rr + 2) % 16];
                2'd2: for (rr = 0; rr < 16; rr = rr + 1) A[5][rr] <= A[5][(rr + 1) % 16];
                2'd3: for (rr = 0; rr < 16; rr = rr + 1) A[5][rr] <= newv[4][(rr + 12) % 16];
            endcase
        end
    end

    always @(posedge clk) begin
        if (state == S_RUN) begin
            case (layer)
                2'd0: for (rr = 0; rr < 16; rr = rr + 1) B[5][rr] <= newv[4][(rr + 1) % 16];
                2'd1: for (rr = 0; rr < 16; rr = rr + 1) B[5][rr] <= newv[4][(rr + 2) % 16];
                2'd2: for (rr = 0; rr < 16; rr = rr + 1) B[5][rr] <= newv[4][(rr + 1) % 16];
                default: ;
            endcase
        end
    end

    // ---- column 6 : inject 6, read 7, edge = 5, present every layer
    always @(posedge clk) begin
        if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd6)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[6][rr] <= (rr == 6) ? inj_val : A[6][(rr + 1) % 16];
        end else if ((state == S_OUT) && (io_cnt[6:4] == 3'd6)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[6][rr] <= A[6][(rr + 1) % 16];
        end else if (commit) begin
            case (layer)
                2'd0: for (rr = 0; rr < 16; rr = rr + 1) A[6][rr] <= A[6][(rr + 3) % 16];
                2'd1: for (rr = 0; rr < 16; rr = rr + 1) A[6][rr] <= A[6][(rr + 14) % 16];
                2'd2: for (rr = 0; rr < 16; rr = rr + 1) A[6][rr] <= A[6][(rr + 14) % 16];
                2'd3: for (rr = 0; rr < 16; rr = rr + 1) A[6][rr] <= newv[5][(rr + 1) % 16];
            endcase
        end
    end

    always @(posedge clk) begin
        if (state == S_RUN) begin
            case (layer)
                2'd0: for (rr = 0; rr < 16; rr = rr + 1) B[6][rr] <= newv[5][(rr + 3) % 16];
                2'd1: for (rr = 0; rr < 16; rr = rr + 1) B[6][rr] <= newv[5][(rr + 14) % 16];
                2'd2: for (rr = 0; rr < 16; rr = rr + 1) B[6][rr] <= newv[5][(rr + 14) % 16];
                default: ;
            endcase
        end
    end

    // ---- column 7 : inject 12, read 13, edge = 6, present every layer
    always @(posedge clk) begin
        if ((state == S_LOAD) && in_data_valid && (io_cnt[6:4] == 3'd7)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[7][rr] <= (rr == 12) ? inj_val : A[7][(rr + 1) % 16];
        end else if ((state == S_OUT) && (io_cnt[6:4] == 3'd7)) begin
            for (rr = 0; rr < 16; rr = rr + 1)
                A[7][rr] <= A[7][(rr + 1) % 16];
        end else if (commit) begin
            case (layer)
                2'd0: for (rr = 0; rr < 16; rr = rr + 1) A[7][rr] <= A[7][(rr + 6) % 16];
                2'd2: for (rr = 0; rr < 16; rr = rr + 1) A[7][rr] <= A[7][(rr + 10) % 16];
                2'd3: for (rr = 0; rr < 16; rr = rr + 1) A[7][rr] <= newv[6][rr];
                default: ; // layer1 self-rotate delta 0 (hold)
            endcase
        end
    end

    always @(posedge clk) begin
        if (state == S_RUN) begin
            case (layer)
                2'd0: for (rr = 0; rr < 16; rr = rr + 1) B[7][rr] <= newv[6][(rr + 6) % 16];
                2'd1: for (rr = 0; rr < 16; rr = rr + 1) B[7][rr] <= newv[6][rr];
                2'd2: for (rr = 0; rr < 16; rr = rr + 1) B[7][rr] <= newv[6][(rr + 10) % 16];
                default: ;
            endcase
        end
    end

    //============================================================
    // Syndrome network: fixed XOR wiring on A's sign bits.
    // rot_amt[n][c] = (BG[n][c] - tf[c]) mod 16, tf = [5,14,10,2,13,12,9,3].
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
    // pr = [11,2,6,14,3,4,7,13] for columns 0..7.
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
