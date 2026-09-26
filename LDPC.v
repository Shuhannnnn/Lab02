module CNU_lane (
    input  [5:0] q_sm0,
    input  [5:0] q_sm1,
    input  [5:0] q_sm2,
    input  [5:0] q_sm3,
    input  [5:0] q_sm4,
    input  [5:0] q_sm5,
    input  [5:0] q_sm6,

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
        mag0 = q_sm0[4:0];
        mag1 = q_sm1[4:0];
        mag2 = q_sm2[4:0];
        mag3 = q_sm3[4:0];
        mag4 = q_sm4[4:0];
        mag5 = q_sm5[4:0];
        mag6 = q_sm6[4:0];

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

        q_neg = {q_sm6[5], q_sm5[5], q_sm4[5], q_sm3[5],
                 q_sm2[5], q_sm1[5], q_sm0[5]};
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

    localparam [1:0] S_IDLE   = 2'd0;
    localparam [1:0] S_READ   = 2'd1;
    localparam [1:0] S_LDPC   = 2'd2;
    localparam [1:0] S_OUTPUT = 2'd3;

    localparam integer CNU_COUNT  = 8;
    localparam integer EDGE_COUNT = 7;

    reg [1:0] state;
    reg [1:0] next_state_c;
    reg       mode_reg;

    reg [6:0] read_idx;
    reg [6:0] output_idx;
    reg [3:0] iteration_cnt;
    reg [2:0] cnu8_cnt;

    wire [1:0] layer;
    wire       batch;
    assign layer = cnu8_cnt[2:1];
    assign batch = cnu8_cnt[0];

    reg signed [7:0] Lja [0:7][0:15];
    reg signed [7:0] Ljb [0:7][0:15];
    reg       [19:0] c2v_queue [0:7][0:7];

    reg        [15:0] pair_vector [0:6][0:7];

    wire syndrome_nonzero;
    wire check_cycle;
    wire at_limit;
    wire continue_iter;
    wire stop_now;
    wire ldpc_commit;
    wire zero_old;
    wire current_bank_sel;
    wire late_input_write;

    reg final_bank_sel;
    reg warn_reg;

    assign check_cycle = (state == S_LDPC) &&
                         (cnu8_cnt == 3'd0) &&
                         (iteration_cnt != 4'd0);
    assign at_limit = (iteration_cnt == 4'd8);
    assign continue_iter = check_cycle && syndrome_nonzero && !at_limit;
    assign stop_now = check_cycle && (!syndrome_nonzero || at_limit);
    assign ldpc_commit = (state == S_LDPC) &&
                         (!check_cycle || continue_iter);
    assign zero_old = (state == S_LDPC) && (iteration_cnt == 4'd0);

    // Mode 0 uses iteration-level ping-pong. Mode 1 always uses Lja.
    assign current_bank_sel = mode_reg ? 1'b0 : iteration_cnt[0];

    // Input[127] is accepted while iteration 1, layer 0, batch 0 is computed.
    assign late_input_write = (state == S_LDPC) &&
                              (iteration_cnt == 4'd0) &&
                              (cnu8_cnt == 3'd0) &&
                              in_data_valid && (read_idx == 7'd127);

    function [15:0] rotate_read16;
        input [15:0] value;
        input [3:0]  shift;
        reg [31:0] doubled;
        begin
            doubled = {value, value};
            rotate_read16 = doubled >> shift;
        end
    endfunction

    //============================================================
    // Fixed QC-H gather routing
    //============================================================

    integer route_lane_i;
    integer route_edge_i;
    always @(*) begin : PAIR_VECTOR
        for (route_edge_i = 0; route_edge_i < EDGE_COUNT;
             route_edge_i = route_edge_i + 1) begin
            for (route_lane_i = 0; route_lane_i < CNU_COUNT;
                 route_lane_i = route_lane_i + 1) begin
                pair_vector[route_edge_i][route_lane_i] = 16'd0;
            end
        end

        case (layer)
            2'd0: begin
                for (route_lane_i = 0; route_lane_i < CNU_COUNT;
                     route_lane_i = route_lane_i + 1) begin
                    pair_vector[0][route_lane_i] = batch ?
                        {Ljb[1][(route_lane_i+6 )&15], Lja[1][(route_lane_i+6 )&15]} :
                        {Ljb[1][(route_lane_i+14)&15], Lja[1][(route_lane_i+14)&15]};
                    pair_vector[1][route_lane_i] = batch ?
                        {Ljb[2][(route_lane_i+2 )&15], Lja[2][(route_lane_i+2 )&15]} :
                        {Ljb[2][(route_lane_i+10)&15], Lja[2][(route_lane_i+10)&15]};
                    pair_vector[2][route_lane_i] = batch ?
                        {Ljb[3][(route_lane_i+10)&15], Lja[3][(route_lane_i+10)&15]} :
                        {Ljb[3][(route_lane_i+2 )&15], Lja[3][(route_lane_i+2 )&15]};
                    pair_vector[3][route_lane_i] = batch ?
                        {Ljb[4][(route_lane_i+5 )&15], Lja[4][(route_lane_i+5 )&15]} :
                        {Ljb[4][(route_lane_i+13)&15], Lja[4][(route_lane_i+13)&15]};
                    pair_vector[4][route_lane_i] = batch ?
                        {Ljb[5][(route_lane_i+4 )&15], Lja[5][(route_lane_i+4 )&15]} :
                        {Ljb[5][(route_lane_i+12)&15], Lja[5][(route_lane_i+12)&15]};
                    pair_vector[5][route_lane_i] = batch ?
                        {Ljb[6][(route_lane_i+1 )&15], Lja[6][(route_lane_i+1 )&15]} :
                        {Ljb[6][(route_lane_i+9 )&15], Lja[6][(route_lane_i+9 )&15]};
                    pair_vector[6][route_lane_i] = batch ?
                        {Ljb[7][(route_lane_i+11)&15], Lja[7][(route_lane_i+11)&15]} :
                        {Ljb[7][(route_lane_i+3 )&15], Lja[7][(route_lane_i+3 )&15]};
                end
            end

            2'd1: begin
                for (route_lane_i = 0; route_lane_i < CNU_COUNT;
                     route_lane_i = route_lane_i + 1) begin
                    pair_vector[0][route_lane_i] = batch ?
                        {Ljb[0][(route_lane_i+13)&15], Lja[0][(route_lane_i+13)&15]} :
                        {Ljb[0][(route_lane_i+5 )&15], Lja[0][(route_lane_i+5 )&15]};
                    pair_vector[1][route_lane_i] = batch ?
                        {Ljb[2][(route_lane_i+6 )&15], Lja[2][(route_lane_i+6 )&15]} :
                        {Ljb[2][(route_lane_i+14)&15], Lja[2][(route_lane_i+14)&15]};
                    pair_vector[2][route_lane_i] = batch ?
                        {Ljb[3][(route_lane_i+2 )&15], Lja[3][(route_lane_i+2 )&15]} :
                        {Ljb[3][(route_lane_i+10)&15], Lja[3][(route_lane_i+10)&15]};
                    pair_vector[3][route_lane_i] = batch ?
                        {Ljb[4][(route_lane_i+10)&15], Lja[4][(route_lane_i+10)&15]} :
                        {Ljb[4][(route_lane_i+2 )&15], Lja[4][(route_lane_i+2 )&15]};
                    pair_vector[4][route_lane_i] = batch ?
                        {Ljb[5][(route_lane_i+5 )&15], Lja[5][(route_lane_i+5 )&15]} :
                        {Ljb[5][(route_lane_i+13)&15], Lja[5][(route_lane_i+13)&15]};
                    pair_vector[5][route_lane_i] = batch ?
                        {Ljb[6][(route_lane_i+4 )&15], Lja[6][(route_lane_i+4 )&15]} :
                        {Ljb[6][(route_lane_i+12)&15], Lja[6][(route_lane_i+12)&15]};
                    pair_vector[6][route_lane_i] = batch ?
                        {Ljb[7][(route_lane_i+1 )&15], Lja[7][(route_lane_i+1 )&15]} :
                        {Ljb[7][(route_lane_i+9 )&15], Lja[7][(route_lane_i+9 )&15]};
                end
            end

            2'd2: begin
                for (route_lane_i = 0; route_lane_i < CNU_COUNT;
                     route_lane_i = route_lane_i + 1) begin
                    pair_vector[0][route_lane_i] = batch ?
                        {Ljb[0][(route_lane_i+8)&15], Lja[0][(route_lane_i+8)&15]} :
                        {Ljb[0][route_lane_i],        Lja[0][route_lane_i]};
                    pair_vector[1][route_lane_i] = batch ?
                        {Ljb[1][(route_lane_i+13)&15], Lja[1][(route_lane_i+13)&15]} :
                        {Ljb[1][(route_lane_i+5 )&15], Lja[1][(route_lane_i+5 )&15]};
                    pair_vector[2][route_lane_i] = batch ?
                        {Ljb[3][(route_lane_i+6 )&15], Lja[3][(route_lane_i+6 )&15]} :
                        {Ljb[3][(route_lane_i+14)&15], Lja[3][(route_lane_i+14)&15]};
                    pair_vector[3][route_lane_i] = batch ?
                        {Ljb[4][(route_lane_i+2 )&15], Lja[4][(route_lane_i+2 )&15]} :
                        {Ljb[4][(route_lane_i+10)&15], Lja[4][(route_lane_i+10)&15]};
                    pair_vector[4][route_lane_i] = batch ?
                        {Ljb[5][(route_lane_i+10)&15], Lja[5][(route_lane_i+10)&15]} :
                        {Ljb[5][(route_lane_i+2 )&15], Lja[5][(route_lane_i+2 )&15]};
                    pair_vector[5][route_lane_i] = batch ?
                        {Ljb[6][(route_lane_i+5 )&15], Lja[6][(route_lane_i+5 )&15]} :
                        {Ljb[6][(route_lane_i+13)&15], Lja[6][(route_lane_i+13)&15]};
                    pair_vector[6][route_lane_i] = batch ?
                        {Ljb[7][(route_lane_i+4 )&15], Lja[7][(route_lane_i+4 )&15]} :
                        {Ljb[7][(route_lane_i+12)&15], Lja[7][(route_lane_i+12)&15]};
                end
            end

            2'd3: begin
                for (route_lane_i = 0; route_lane_i < CNU_COUNT;
                     route_lane_i = route_lane_i + 1) begin
                    pair_vector[0][route_lane_i] = batch ?
                        {Ljb[0][(route_lane_i+15)&15], Lja[0][(route_lane_i+15)&15]} :
                        {Ljb[0][(route_lane_i+7 )&15], Lja[0][(route_lane_i+7 )&15]};
                    pair_vector[1][route_lane_i] = batch ?
                        {Ljb[1][(route_lane_i+8)&15], Lja[1][(route_lane_i+8)&15]} :
                        {Ljb[1][route_lane_i],       Lja[1][route_lane_i]};
                    pair_vector[2][route_lane_i] = batch ?
                        {Ljb[2][(route_lane_i+13)&15], Lja[2][(route_lane_i+13)&15]} :
                        {Ljb[2][(route_lane_i+5 )&15], Lja[2][(route_lane_i+5 )&15]};
                    pair_vector[3][route_lane_i] = batch ?
                        {Ljb[4][(route_lane_i+6 )&15], Lja[4][(route_lane_i+6 )&15]} :
                        {Ljb[4][(route_lane_i+14)&15], Lja[4][(route_lane_i+14)&15]};
                    pair_vector[4][route_lane_i] = batch ?
                        {Ljb[5][(route_lane_i+2 )&15], Lja[5][(route_lane_i+2 )&15]} :
                        {Ljb[5][(route_lane_i+10)&15], Lja[5][(route_lane_i+10)&15]};
                    pair_vector[5][route_lane_i] = batch ?
                        {Ljb[6][(route_lane_i+10)&15], Lja[6][(route_lane_i+10)&15]} :
                        {Ljb[6][(route_lane_i+2 )&15], Lja[6][(route_lane_i+2 )&15]};
                    pair_vector[6][route_lane_i] = batch ?
                        {Ljb[7][(route_lane_i+5 )&15], Lja[7][(route_lane_i+5 )&15]} :
                        {Ljb[7][(route_lane_i+13)&15], Lja[7][(route_lane_i+13)&15]};
                end
            end

            default: begin
            end
        endcase
    end

    //============================================================
    // FSM and frame control
    //============================================================

    always @(*) begin : STATE_NEXT
        next_state_c = state;
        case (state)
            S_IDLE: begin
                if (in_mode_valid && !out_valid)
                    next_state_c = S_READ;
            end
            S_READ: begin
                if (in_data_valid && (read_idx == 7'd126))
                    next_state_c = S_LDPC;
            end
            S_LDPC: begin
                if (stop_now)
                    next_state_c = S_OUTPUT;
            end
            S_OUTPUT: begin
                if (output_idx == 7'd127)
                    next_state_c = S_IDLE;
            end
            default: begin
                next_state_c = S_IDLE;
            end
        endcase
    end

    always @(posedge clk or negedge rst_n) begin : FSM_STATE
        if (!rst_n)
            state <= S_IDLE;
        else
            state <= next_state_c;
    end

    always @(posedge clk or negedge rst_n) begin : MODE_CONTROL
        if (!rst_n)
            mode_reg <= 1'b0;
        else if (in_mode_valid && (state == S_IDLE) && !out_valid)
            mode_reg <= in_mode;
    end

    always @(posedge clk or negedge rst_n) begin : READ_INDEX
        if (!rst_n) begin
            read_idx <= 7'd0;
        end else if (in_mode_valid && (state == S_IDLE) && !out_valid) begin
            read_idx <= 7'd0;
        end else if ((state == S_READ) && in_data_valid) begin
            read_idx <= read_idx + 7'd1;
        end else if (late_input_write) begin
            read_idx <= 7'd0;
        end
    end

    always @(posedge clk or negedge rst_n) begin : LDPC_COUNTER
        if (!rst_n) begin
            iteration_cnt <= 4'd0;
            cnu8_cnt <= 3'd0;
        end else if (in_mode_valid && (state == S_IDLE) && !out_valid) begin
            iteration_cnt <= 4'd0;
            cnu8_cnt <= 3'd0;
        end else if (state == S_LDPC) begin
            if (stop_now) begin
                cnu8_cnt <= 3'd0;
            end else if (ldpc_commit) begin
                if (cnu8_cnt == 3'd7) begin
                    cnu8_cnt <= 3'd0;
                    iteration_cnt <= iteration_cnt + 4'd1;
                end else begin
                    cnu8_cnt <= cnu8_cnt + 3'd1;
                end
            end
        end
    end

    //============================================================
    // Old C2V reconstruction and edge arithmetic
    //============================================================

    reg        [4:0] min1_old    [0:7];
    reg        [4:0] min2_old    [0:7];
    reg        [2:0] min1_idx    [0:7];
    reg        [6:0] c2v_neg     [0:7];
    reg        [4:0] r_old_mag   [0:7][0:6];
    reg        [7:0] r_sub_operand [0:7][0:6];
    reg              r_sub_carry   [0:7][0:6];
    reg        [8:0] a_base_sum  [0:7][0:6];
    reg        [8:0] b_base_sum  [0:7][0:6];
    reg signed [7:0] a_base      [0:7][0:6];
    reg signed [7:0] b_base      [0:7][0:6];
    reg signed [7:0] q_base      [0:7][0:6];
    reg        [5:0] q_sm        [0:7][0:6];
    reg signed [7:0] update_base [0:7][0:6];
    reg              use_q_base  [0:6];
    reg              update_bank_sel [0:6];

    reg [6:0] first_touch_edge;
    always @(*) begin : FIRST_TOUCH
        case (layer)
            2'd0: first_touch_edge = 7'b1111111;
            2'd1: first_touch_edge = 7'b0000001;
            default: first_touch_edge = 7'b0000000;
        endcase
    end

    integer cnu_i;
    integer edge_i;
    always @(*) begin : CNU_INPUT
        for (edge_i = 0; edge_i < EDGE_COUNT; edge_i = edge_i + 1) begin
            use_q_base[edge_i] = mode_reg || first_touch_edge[edge_i];
            update_bank_sel[edge_i] = current_bank_sel ^
                                      ~use_q_base[edge_i];
        end

        for (cnu_i = 0; cnu_i < CNU_COUNT; cnu_i = cnu_i + 1) begin
            if (old_c2v_valid) begin
                min1_old[cnu_i] = c2v_queue[0][cnu_i][19:15];
                min2_old[cnu_i] = c2v_queue[0][cnu_i][14:10];
                min1_idx[cnu_i] = c2v_queue[0][cnu_i][9:7];
                c2v_neg[cnu_i]  = c2v_queue[0][cnu_i][6:0];
            end else begin
                min1_old[cnu_i] = 5'd0;
                min2_old[cnu_i] = 5'd0;
                min1_idx[cnu_i] = 3'd0;
                c2v_neg[cnu_i]  = 7'd0;
            end

<<<<<<< HEAD
        if (state == S_LDPC) begin
            for (cnu_i = 0; cnu_i < CNU_COUNT; cnu_i = cnu_i + 1) begin
                if (!zero_old) begin
                    min1_old[cnu_i] = c2v[layer][batch][cnu_i][19:15];
                    min2_old[cnu_i] = c2v[layer][batch][cnu_i][14:10];
                    min1_idx[cnu_i] = c2v[layer][batch][cnu_i][9:7];
                    c2v_neg[cnu_i]  = c2v[layer][batch][cnu_i][6:0];
                end

                for (edge_i = 0; edge_i < EDGE_COUNT; edge_i = edge_i + 1) begin
                    if (zero_old) begin
                        r_old_mag[cnu_i][edge_i] = 5'd0;
                        r_old[cnu_i][edge_i] = 6'sd0;
                    end else begin
                        if (min1_idx[cnu_i] == edge_i)
                            r_old_mag[cnu_i][edge_i] = min2_old[cnu_i];
                        else
                            r_old_mag[cnu_i][edge_i] = min1_old[cnu_i];

                        if (c2v_neg[cnu_i][edge_i])
                            r_old[cnu_i][edge_i] =
                                -$signed({1'b0, r_old_mag[cnu_i][edge_i]});
                        else
                            r_old[cnu_i][edge_i] =
                                 $signed({1'b0, r_old_mag[cnu_i][edge_i]});
                    end

                    r_old_ext[cnu_i][edge_i] =
                        {{2{r_old[cnu_i][edge_i][5]}}, r_old[cnu_i][edge_i]};
                    q_base[cnu_i][edge_i] =
                        src_vector[edge_i][cnu_i] - r_old_ext[cnu_i][edge_i];
                    dst_base[cnu_i][edge_i] =
                        dst_vector[edge_i][cnu_i] - r_old_ext[cnu_i][edge_i];
                    q_msg[cnu_i][edge_i] = clip6(q_base[cnu_i][edge_i]);

                    if (mode_reg || first_touch_edge[edge_i])
                        update_base[cnu_i][edge_i] = q_base[cnu_i][edge_i];
                    else
                        update_base[cnu_i][edge_i] = dst_base[cnu_i][edge_i];
=======
            for (edge_i = 0; edge_i < EDGE_COUNT; edge_i = edge_i + 1) begin
                if (min1_idx[cnu_i] == edge_i)
                    r_old_mag[cnu_i][edge_i] = min2_old[cnu_i];
                else
                    r_old_mag[cnu_i][edge_i] = min1_old[cnu_i];

                // Subtract R_old from both physical banks before selecting
                // either the source value or the accumulation destination.
                r_sub_carry[cnu_i][edge_i] =
                    ~c2v_neg[cnu_i][edge_i];
                r_sub_operand[cnu_i][edge_i] =
                    {3'b000, r_old_mag[cnu_i][edge_i]} ^
                    {8{r_sub_carry[cnu_i][edge_i]}};

                a_base_sum[cnu_i][edge_i] =
                    {1'b0, pair_vector[edge_i][cnu_i][7:0]} +
                    {1'b0, r_sub_operand[cnu_i][edge_i]} +
                    r_sub_carry[cnu_i][edge_i];
                b_base_sum[cnu_i][edge_i] =
                    {1'b0, pair_vector[edge_i][cnu_i][15:8]} +
                    {1'b0, r_sub_operand[cnu_i][edge_i]} +
                    r_sub_carry[cnu_i][edge_i];
                a_base[cnu_i][edge_i] =
                    $signed(a_base_sum[cnu_i][edge_i][7:0]);
                b_base[cnu_i][edge_i] =
                    $signed(b_base_sum[cnu_i][edge_i][7:0]);

                if (current_bank_sel)
                    q_base[cnu_i][edge_i] = b_base[cnu_i][edge_i];
                else
                    q_base[cnu_i][edge_i] = a_base[cnu_i][edge_i];

                if (update_bank_sel[edge_i])
                    update_base[cnu_i][edge_i] = b_base[cnu_i][edge_i];
                else
                    update_base[cnu_i][edge_i] = a_base[cnu_i][edge_i];

                // clip6(q_base) followed by abs6 is represented directly as
                // one sign bit and one saturated five-bit magnitude.
                q_sm[cnu_i][edge_i][5] = q_base[cnu_i][edge_i][7];
                if ((q_base[cnu_i][edge_i][6] ^
                     q_base[cnu_i][edge_i][7]) ||
                    (q_base[cnu_i][edge_i][5] ^
                     q_base[cnu_i][edge_i][7]) ||
                    (q_base[cnu_i][edge_i][7] &&
                     (q_base[cnu_i][edge_i][4:0] == 5'd0))) begin
                    q_sm[cnu_i][edge_i][4:0] = 5'd31;
                end else if (q_base[cnu_i][edge_i][7]) begin
                    q_sm[cnu_i][edge_i][4:0] =
                        (~q_base[cnu_i][edge_i][4:0]) + 5'd1;
                end else begin
                    q_sm[cnu_i][edge_i][4:0] =
                        q_base[cnu_i][edge_i][4:0];
>>>>>>> 8CNU_03fail
                end
            end
        end
    end

    wire        [19:0] c2v_new_lane    [0:7];
    wire        [41:0] r_new_flat_lane [0:7];
    wire signed  [5:0] r_new            [0:7][0:6];
    wire signed  [7:0] r_new_ext        [0:7][0:6];
    wire signed  [7:0] new_lj           [0:7][0:6];

    genvar cnu_g;
    genvar edge_g;
    generate
        for (cnu_g = 0; cnu_g < CNU_COUNT; cnu_g = cnu_g + 1) begin : GEN_CNU
            CNU_lane u_CNU_lane (
                .q_sm0(q_sm[cnu_g][0]),
                .q_sm1(q_sm[cnu_g][1]),
                .q_sm2(q_sm[cnu_g][2]),
                .q_sm3(q_sm[cnu_g][3]),
                .q_sm4(q_sm[cnu_g][4]),
                .q_sm5(q_sm[cnu_g][5]),
                .q_sm6(q_sm[cnu_g][6]),
                .c2v_new(c2v_new_lane[cnu_g]),
                .r_new_flat(r_new_flat_lane[cnu_g])
            );

            for (edge_g = 0; edge_g < EDGE_COUNT; edge_g = edge_g + 1) begin : GEN_NEW_LJ
                assign r_new[cnu_g][edge_g] =
                    r_new_flat_lane[cnu_g][edge_g*6 +: 6];
                assign r_new_ext[cnu_g][edge_g] =
                    {{2{r_new[cnu_g][edge_g][5]}}, r_new[cnu_g][edge_g]};
                assign new_lj[cnu_g][edge_g] =
                    $signed(update_base[cnu_g][edge_g]) +
                    $signed(r_new_ext[cnu_g][edge_g]);
            end
        end
    endgenerate

    //============================================================
    // Posterior and descriptor writeback
    //============================================================

    function [2:0] h_edge_for_col;
        input [1:0] layer_i;
        input [2:0] col_i;
        begin
            if (col_i == {1'b0, layer_i})
                h_edge_for_col = 3'd0;
            else if (col_i > {1'b0, layer_i})
                h_edge_for_col = col_i - 3'd1;
            else
                h_edge_for_col = col_i;
        end
    endfunction

    function [3:0] h_shift_for_col;
        input [1:0] layer_i;
        input [2:0] col_i;
        begin
            h_shift_for_col = 4'd0;
            case ({layer_i, col_i})
                {2'd0, 3'd1}: h_shift_for_col = 4'd14;
                {2'd0, 3'd2}: h_shift_for_col = 4'd10;
                {2'd0, 3'd3}: h_shift_for_col = 4'd2;
                {2'd0, 3'd4}: h_shift_for_col = 4'd13;
                {2'd0, 3'd5}: h_shift_for_col = 4'd12;
                {2'd0, 3'd6}: h_shift_for_col = 4'd9;
                {2'd0, 3'd7}: h_shift_for_col = 4'd3;

                {2'd1, 3'd0}: h_shift_for_col = 4'd5;
                {2'd1, 3'd2}: h_shift_for_col = 4'd14;
                {2'd1, 3'd3}: h_shift_for_col = 4'd10;
                {2'd1, 3'd4}: h_shift_for_col = 4'd2;
                {2'd1, 3'd5}: h_shift_for_col = 4'd13;
                {2'd1, 3'd6}: h_shift_for_col = 4'd12;
                {2'd1, 3'd7}: h_shift_for_col = 4'd9;

                {2'd2, 3'd0}: h_shift_for_col = 4'd0;
                {2'd2, 3'd1}: h_shift_for_col = 4'd5;
                {2'd2, 3'd3}: h_shift_for_col = 4'd14;
                {2'd2, 3'd4}: h_shift_for_col = 4'd10;
                {2'd2, 3'd5}: h_shift_for_col = 4'd2;
                {2'd2, 3'd6}: h_shift_for_col = 4'd13;
                {2'd2, 3'd7}: h_shift_for_col = 4'd12;

                {2'd3, 3'd0}: h_shift_for_col = 4'd7;
                {2'd3, 3'd1}: h_shift_for_col = 4'd0;
                {2'd3, 3'd2}: h_shift_for_col = 4'd5;
                {2'd3, 3'd4}: h_shift_for_col = 4'd14;
                {2'd3, 3'd5}: h_shift_for_col = 4'd10;
                {2'd3, 3'd6}: h_shift_for_col = 4'd2;
                {2'd3, 3'd7}: h_shift_for_col = 4'd13;
                default: h_shift_for_col = 4'd0;
            endcase
        end
    endfunction

    reg  [7:0]  phase_sel;
    reg  [7:0]  read_col_sel;
    reg  [15:0] read_pos_sel;
    wire         write_lja;

    assign write_lja = mode_reg || current_bank_sel;

    always @(*) begin : SHARED_WRITE_DECODE
        phase_sel = 8'd0;
        phase_sel[cnu8_cnt] = 1'b1;
        read_col_sel = 8'd0;
        read_col_sel[read_idx[6:4]] = 1'b1;
        read_pos_sel = 16'd0;
        read_pos_sel[read_idx[3:0]] = 1'b1;
    end

    genvar wb_col_g;
    genvar wb_pos_g;
    generate
        for (wb_col_g = 0; wb_col_g < 8;
             wb_col_g = wb_col_g + 1) begin : GEN_WB_COL
            for (wb_pos_g = 0; wb_pos_g < 16;
                 wb_pos_g = wb_pos_g + 1) begin : GEN_WB_POS
                localparam [2:0] WB_EDGE0 =
                    h_edge_for_col(2'd0, wb_col_g);
                localparam [2:0] WB_EDGE1 =
                    h_edge_for_col(2'd1, wb_col_g);
                localparam [2:0] WB_EDGE2 =
                    h_edge_for_col(2'd2, wb_col_g);
                localparam [2:0] WB_EDGE3 =
                    h_edge_for_col(2'd3, wb_col_g);

                localparam [3:0] WB_SHIFT0 =
                    h_shift_for_col(2'd0, wb_col_g);
                localparam [3:0] WB_SHIFT1 =
                    h_shift_for_col(2'd1, wb_col_g);
                localparam [3:0] WB_SHIFT2 =
                    h_shift_for_col(2'd2, wb_col_g);
                localparam [3:0] WB_SHIFT3 =
                    h_shift_for_col(2'd3, wb_col_g);

                localparam integer WB_ROW0 =
                    (wb_pos_g + 16 - WB_SHIFT0) & 15;
                localparam integer WB_ROW1 =
                    (wb_pos_g + 16 - WB_SHIFT1) & 15;
                localparam integer WB_ROW2 =
                    (wb_pos_g + 16 - WB_SHIFT2) & 15;
                localparam integer WB_ROW3 =
                    (wb_pos_g + 16 - WB_SHIFT3) & 15;

                localparam integer WB_LANE0 = WB_ROW0 & 7;
                localparam integer WB_LANE1 = WB_ROW1 & 7;
                localparam integer WB_LANE2 = WB_ROW2 & 7;
                localparam integer WB_LANE3 = WB_ROW3 & 7;

                localparam integer WB_PHASE0 = (WB_ROW0 >> 3);
                localparam integer WB_PHASE1 = 2 + (WB_ROW1 >> 3);
                localparam integer WB_PHASE2 = 4 + (WB_ROW2 >> 3);
                localparam integer WB_PHASE3 = 6 + (WB_ROW3 >> 3);

                wire wb_sel0;
                wire wb_sel1;
                wire wb_sel2;
                wire wb_sel3;
                wire wb_hit;
                wire signed [7:0] wb_data;

                assign wb_sel0 = (wb_col_g != 0) ?
                                 phase_sel[WB_PHASE0] : 1'b0;
                assign wb_sel1 = (wb_col_g != 1) ?
                                 phase_sel[WB_PHASE1] : 1'b0;
                assign wb_sel2 = (wb_col_g != 2) ?
                                 phase_sel[WB_PHASE2] : 1'b0;
                assign wb_sel3 = (wb_col_g != 3) ?
                                 phase_sel[WB_PHASE3] : 1'b0;
                assign wb_hit = wb_sel0 || wb_sel1 || wb_sel2 || wb_sel3;

                assign wb_data =
                    ({8{wb_sel0}} & new_lj[WB_LANE0][WB_EDGE0]) |
                    ({8{wb_sel1}} & new_lj[WB_LANE1][WB_EDGE1]) |
                    ({8{wb_sel2}} & new_lj[WB_LANE2][WB_EDGE2]) |
                    ({8{wb_sel3}} & new_lj[WB_LANE3][WB_EDGE3]);

                // Each physical Lj word has exactly one sequential owner.
                // Priority matches the original nonblocking assignment order.
                always @(posedge clk) begin
                    if (late_input_write &&
                        (wb_col_g == 7) && (wb_pos_g == 15)) begin
                        Lja[wb_col_g][wb_pos_g] <=
                            {{2{in_data[5]}}, in_data};
                    end else if (ldpc_commit && wb_hit && write_lja) begin
                        Lja[wb_col_g][wb_pos_g] <= wb_data;
                    end else if ((state == S_READ) && in_data_valid &&
                                 read_col_sel[wb_col_g] &&
                                 read_pos_sel[wb_pos_g]) begin
                        Lja[wb_col_g][wb_pos_g] <=
                            {{2{in_data[5]}}, in_data};
                    end

                    if (ldpc_commit && wb_hit && !write_lja)
                        Ljb[wb_col_g][wb_pos_g] <= wb_data;
                end
            end
        end
    endgenerate

    integer c2v_stage_i;
    integer c2v_lane_i;
    always @(posedge clk) begin : C2V_QUEUE_WRITEBACK
        if (ldpc_commit) begin
            for (c2v_stage_i = 0; c2v_stage_i < 7;
                 c2v_stage_i = c2v_stage_i + 1) begin
                for (c2v_lane_i = 0; c2v_lane_i < CNU_COUNT;
                     c2v_lane_i = c2v_lane_i + 1) begin
                    c2v_queue[c2v_stage_i][c2v_lane_i] <=
                        c2v_queue[c2v_stage_i+1][c2v_lane_i];
                end
            end
            for (c2v_lane_i = 0; c2v_lane_i < CNU_COUNT;
                 c2v_lane_i = c2v_lane_i + 1) begin
                c2v_queue[7][c2v_lane_i] <=
                    c2v_new_lane[c2v_lane_i];
            end
        end
    end

    //============================================================
    // Fixed QC syndrome network
    //============================================================

    reg [15:0] hard0;
    reg [15:0] hard1;
    reg [15:0] hard2;
    reg [15:0] hard3;
    reg [15:0] hard4;
    reg [15:0] hard5;
    reg [15:0] hard6;
    reg [15:0] hard7;

    integer hard_i;
    always @(*) begin : HARD_DECISION
        for (hard_i = 0; hard_i < 16; hard_i = hard_i + 1) begin
            if (current_bank_sel) begin
                hard0[hard_i] = Ljb[0][hard_i][7];
                hard1[hard_i] = Ljb[1][hard_i][7];
                hard2[hard_i] = Ljb[2][hard_i][7];
                hard3[hard_i] = Ljb[3][hard_i][7];
                hard4[hard_i] = Ljb[4][hard_i][7];
                hard5[hard_i] = Ljb[5][hard_i][7];
                hard6[hard_i] = Ljb[6][hard_i][7];
                hard7[hard_i] = Ljb[7][hard_i][7];
            end else begin
                hard0[hard_i] = Lja[0][hard_i][7];
                hard1[hard_i] = Lja[1][hard_i][7];
                hard2[hard_i] = Lja[2][hard_i][7];
                hard3[hard_i] = Lja[3][hard_i][7];
                hard4[hard_i] = Lja[4][hard_i][7];
                hard5[hard_i] = Lja[5][hard_i][7];
                hard6[hard_i] = Lja[6][hard_i][7];
                hard7[hard_i] = Lja[7][hard_i][7];
            end
        end
    end

    wire [15:0] syn0_pair0;
    wire [15:0] syn0_pair1;
    wire [15:0] syn0_pair2;
    wire [15:0] syn1_pair0;
    wire [15:0] syn1_pair1;
    wire [15:0] syn1_pair2;
    wire [15:0] syn2_pair0;
    wire [15:0] syn2_pair1;
    wire [15:0] syn2_pair2;
    wire [15:0] syn3_pair0;
    wire [15:0] syn3_pair1;
    wire [15:0] syn3_pair2;
    wire [15:0] syn_vec0;
    wire [15:0] syn_vec1;
    wire [15:0] syn_vec2;
    wire [15:0] syn_vec3;

    assign syn0_pair0 = rotate_read16(hard1, 4'd14) ^ rotate_read16(hard2, 4'd10);
    assign syn0_pair1 = rotate_read16(hard3, 4'd2 ) ^ rotate_read16(hard4, 4'd13);
    assign syn0_pair2 = rotate_read16(hard5, 4'd12) ^ rotate_read16(hard6, 4'd9 );
    assign syn_vec0 = (syn0_pair0 ^ syn0_pair1) ^
                      (syn0_pair2 ^ rotate_read16(hard7, 4'd3));

    assign syn1_pair0 = rotate_read16(hard0, 4'd5 ) ^ rotate_read16(hard2, 4'd14);
    assign syn1_pair1 = rotate_read16(hard3, 4'd10) ^ rotate_read16(hard4, 4'd2 );
    assign syn1_pair2 = rotate_read16(hard5, 4'd13) ^ rotate_read16(hard6, 4'd12);
    assign syn_vec1 = (syn1_pair0 ^ syn1_pair1) ^
                      (syn1_pair2 ^ rotate_read16(hard7, 4'd9));

    assign syn2_pair0 = hard0 ^ rotate_read16(hard1, 4'd5);
    assign syn2_pair1 = rotate_read16(hard3, 4'd14) ^ rotate_read16(hard4, 4'd10);
    assign syn2_pair2 = rotate_read16(hard5, 4'd2 ) ^ rotate_read16(hard6, 4'd13);
    assign syn_vec2 = (syn2_pair0 ^ syn2_pair1) ^
                      (syn2_pair2 ^ rotate_read16(hard7, 4'd12));

    assign syn3_pair0 = rotate_read16(hard0, 4'd7) ^ hard1;
    assign syn3_pair1 = rotate_read16(hard2, 4'd5 ) ^ rotate_read16(hard4, 4'd14);
    assign syn3_pair2 = rotate_read16(hard5, 4'd10) ^ rotate_read16(hard6, 4'd2 );
    assign syn_vec3 = (syn3_pair0 ^ syn3_pair1) ^
                      (syn3_pair2 ^ rotate_read16(hard7, 4'd13));

    assign syndrome_nonzero = |(syn_vec0 | syn_vec1 | syn_vec2 | syn_vec3);

    //============================================================
    // Registered output
    //============================================================

    always @(posedge clk or negedge rst_n) begin : OUTPUT_CONTROL
        if (!rst_n) begin
            out_valid <= 1'b0;
            out_data <= 8'd0;
            out_warn <= 1'b0;
            output_idx <= 7'd0;
            final_bank_sel <= 1'b0;
            warn_reg <= 1'b0;
        end else if (stop_now) begin
            final_bank_sel <= current_bank_sel;
            warn_reg <= at_limit && syndrome_nonzero;
            out_valid <= 1'b1;
            out_warn <= at_limit && syndrome_nonzero;
            if (current_bank_sel)
                out_data <= Ljb[0][0];
            else
                out_data <= Lja[0][0];
            output_idx <= 7'd1;
        end else if (state == S_OUTPUT) begin
            out_valid <= 1'b1;
            out_warn <= warn_reg;
            if (final_bank_sel)
                out_data <= Ljb[output_idx[6:4]][output_idx[3:0]];
            else
                out_data <= Lja[output_idx[6:4]][output_idx[3:0]];

            if (output_idx == 7'd127)
                output_idx <= 7'd0;
            else
                output_idx <= output_idx + 7'd1;
        end else begin
            out_valid <= 1'b0;
            out_data <= 8'd0;
            out_warn <= 1'b0;
            output_idx <= 7'd0;
        end
    end

endmodule
