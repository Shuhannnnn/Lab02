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
    reg       [19:0] c2v [0:3][0:1][0:7];

    reg        [15:0] pair_vector [0:6][0:7];
    reg signed  [7:0] src_vector  [0:6][0:7];
    reg signed  [7:0] dst_vector  [0:6][0:7];

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

    integer src_dst_edge_i;
    integer src_dst_lane_i;
    always @(*) begin : SRC_DST_ROUTING
        for (src_dst_edge_i = 0; src_dst_edge_i < EDGE_COUNT;
             src_dst_edge_i = src_dst_edge_i + 1) begin
            for (src_dst_lane_i = 0; src_dst_lane_i < CNU_COUNT;
                 src_dst_lane_i = src_dst_lane_i + 1) begin
                if (current_bank_sel == 1'b0) begin
                    src_vector[src_dst_edge_i][src_dst_lane_i] =
                        $signed(pair_vector[src_dst_edge_i][src_dst_lane_i][7:0]);
                    dst_vector[src_dst_edge_i][src_dst_lane_i] =
                        $signed(pair_vector[src_dst_edge_i][src_dst_lane_i][15:8]);
                end else begin
                    src_vector[src_dst_edge_i][src_dst_lane_i] =
                        $signed(pair_vector[src_dst_edge_i][src_dst_lane_i][15:8]);
                    dst_vector[src_dst_edge_i][src_dst_lane_i] =
                        $signed(pair_vector[src_dst_edge_i][src_dst_lane_i][7:0]);
                end
            end
        end
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

    reg signed [5:0] r_old       [0:7][0:6];
    reg        [4:0] min1_old    [0:7];
    reg        [4:0] min2_old    [0:7];
    reg        [2:0] min1_idx    [0:7];
    reg        [6:0] c2v_neg     [0:7];
    reg        [4:0] r_old_mag   [0:7][0:6];
    reg signed [7:0] r_old_ext   [0:7][0:6];
    reg signed [7:0] q_base      [0:7][0:6];
    reg signed [7:0] dst_base    [0:7][0:6];
    reg signed [5:0] q_msg       [0:7][0:6];
    reg signed [7:0] update_base [0:7][0:6];

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
        for (cnu_i = 0; cnu_i < CNU_COUNT; cnu_i = cnu_i + 1) begin
            min1_old[cnu_i] = 5'd0;
            min2_old[cnu_i] = 5'd0;
            min1_idx[cnu_i] = 3'd0;
            c2v_neg[cnu_i]  = 7'd0;
            for (edge_i = 0; edge_i < EDGE_COUNT; edge_i = edge_i + 1) begin
                r_old_mag[cnu_i][edge_i] = 5'd0;
                r_old[cnu_i][edge_i] = 6'sd0;
                r_old_ext[cnu_i][edge_i] = 8'sd0;
                q_base[cnu_i][edge_i] = 8'sd0;
                dst_base[cnu_i][edge_i] = 8'sd0;
                q_msg[cnu_i][edge_i] = 6'sd0;
                update_base[cnu_i][edge_i] = 8'sd0;
            end
        end

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
                .q0(q_msg[cnu_g][0]),
                .q1(q_msg[cnu_g][1]),
                .q2(q_msg[cnu_g][2]),
                .q3(q_msg[cnu_g][3]),
                .q4(q_msg[cnu_g][4]),
                .q5(q_msg[cnu_g][5]),
                .q6(q_msg[cnu_g][6]),
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

    integer lj_i;
    always @(posedge clk) begin : POSTERIOR_WRITEBACK
        if ((state == S_READ) && in_data_valid) begin
            Lja[read_idx[6:4]][read_idx[3:0]] <=
                {{2{in_data[5]}}, in_data};
        end

        if (ldpc_commit) begin
            case (layer)
                2'd0: begin
                    for (lj_i = 0; lj_i < CNU_COUNT; lj_i = lj_i + 1) begin
                        if (mode_reg || current_bank_sel) begin
                            Lja[1][batch ? ((lj_i+6 )&15) : ((lj_i+14)&15)] <= new_lj[lj_i][0];
                            Lja[2][batch ? ((lj_i+2 )&15) : ((lj_i+10)&15)] <= new_lj[lj_i][1];
                            Lja[3][batch ? ((lj_i+10)&15) : ((lj_i+2 )&15)] <= new_lj[lj_i][2];
                            Lja[4][batch ? ((lj_i+5 )&15) : ((lj_i+13)&15)] <= new_lj[lj_i][3];
                            Lja[5][batch ? ((lj_i+4 )&15) : ((lj_i+12)&15)] <= new_lj[lj_i][4];
                            Lja[6][batch ? ((lj_i+1 )&15) : ((lj_i+9 )&15)] <= new_lj[lj_i][5];
                            Lja[7][batch ? ((lj_i+11)&15) : ((lj_i+3 )&15)] <= new_lj[lj_i][6];
                        end else begin
                            Ljb[1][batch ? ((lj_i+6 )&15) : ((lj_i+14)&15)] <= new_lj[lj_i][0];
                            Ljb[2][batch ? ((lj_i+2 )&15) : ((lj_i+10)&15)] <= new_lj[lj_i][1];
                            Ljb[3][batch ? ((lj_i+10)&15) : ((lj_i+2 )&15)] <= new_lj[lj_i][2];
                            Ljb[4][batch ? ((lj_i+5 )&15) : ((lj_i+13)&15)] <= new_lj[lj_i][3];
                            Ljb[5][batch ? ((lj_i+4 )&15) : ((lj_i+12)&15)] <= new_lj[lj_i][4];
                            Ljb[6][batch ? ((lj_i+1 )&15) : ((lj_i+9 )&15)] <= new_lj[lj_i][5];
                            Ljb[7][batch ? ((lj_i+11)&15) : ((lj_i+3 )&15)] <= new_lj[lj_i][6];
                        end
                    end
                end

                2'd1: begin
                    for (lj_i = 0; lj_i < CNU_COUNT; lj_i = lj_i + 1) begin
                        if (mode_reg || current_bank_sel) begin
                            Lja[0][batch ? ((lj_i+13)&15) : ((lj_i+5 )&15)] <= new_lj[lj_i][0];
                            Lja[2][batch ? ((lj_i+6 )&15) : ((lj_i+14)&15)] <= new_lj[lj_i][1];
                            Lja[3][batch ? ((lj_i+2 )&15) : ((lj_i+10)&15)] <= new_lj[lj_i][2];
                            Lja[4][batch ? ((lj_i+10)&15) : ((lj_i+2 )&15)] <= new_lj[lj_i][3];
                            Lja[5][batch ? ((lj_i+5 )&15) : ((lj_i+13)&15)] <= new_lj[lj_i][4];
                            Lja[6][batch ? ((lj_i+4 )&15) : ((lj_i+12)&15)] <= new_lj[lj_i][5];
                            Lja[7][batch ? ((lj_i+1 )&15) : ((lj_i+9 )&15)] <= new_lj[lj_i][6];
                        end else begin
                            Ljb[0][batch ? ((lj_i+13)&15) : ((lj_i+5 )&15)] <= new_lj[lj_i][0];
                            Ljb[2][batch ? ((lj_i+6 )&15) : ((lj_i+14)&15)] <= new_lj[lj_i][1];
                            Ljb[3][batch ? ((lj_i+2 )&15) : ((lj_i+10)&15)] <= new_lj[lj_i][2];
                            Ljb[4][batch ? ((lj_i+10)&15) : ((lj_i+2 )&15)] <= new_lj[lj_i][3];
                            Ljb[5][batch ? ((lj_i+5 )&15) : ((lj_i+13)&15)] <= new_lj[lj_i][4];
                            Ljb[6][batch ? ((lj_i+4 )&15) : ((lj_i+12)&15)] <= new_lj[lj_i][5];
                            Ljb[7][batch ? ((lj_i+1 )&15) : ((lj_i+9 )&15)] <= new_lj[lj_i][6];
                        end
                    end
                end

                2'd2: begin
                    for (lj_i = 0; lj_i < CNU_COUNT; lj_i = lj_i + 1) begin
                        if (mode_reg || current_bank_sel) begin
                            Lja[0][batch ? ((lj_i+8 )&15) : lj_i] <= new_lj[lj_i][0];
                            Lja[1][batch ? ((lj_i+13)&15) : ((lj_i+5 )&15)] <= new_lj[lj_i][1];
                            Lja[3][batch ? ((lj_i+6 )&15) : ((lj_i+14)&15)] <= new_lj[lj_i][2];
                            Lja[4][batch ? ((lj_i+2 )&15) : ((lj_i+10)&15)] <= new_lj[lj_i][3];
                            Lja[5][batch ? ((lj_i+10)&15) : ((lj_i+2 )&15)] <= new_lj[lj_i][4];
                            Lja[6][batch ? ((lj_i+5 )&15) : ((lj_i+13)&15)] <= new_lj[lj_i][5];
                            Lja[7][batch ? ((lj_i+4 )&15) : ((lj_i+12)&15)] <= new_lj[lj_i][6];
                        end else begin
                            Ljb[0][batch ? ((lj_i+8 )&15) : lj_i] <= new_lj[lj_i][0];
                            Ljb[1][batch ? ((lj_i+13)&15) : ((lj_i+5 )&15)] <= new_lj[lj_i][1];
                            Ljb[3][batch ? ((lj_i+6 )&15) : ((lj_i+14)&15)] <= new_lj[lj_i][2];
                            Ljb[4][batch ? ((lj_i+2 )&15) : ((lj_i+10)&15)] <= new_lj[lj_i][3];
                            Ljb[5][batch ? ((lj_i+10)&15) : ((lj_i+2 )&15)] <= new_lj[lj_i][4];
                            Ljb[6][batch ? ((lj_i+5 )&15) : ((lj_i+13)&15)] <= new_lj[lj_i][5];
                            Ljb[7][batch ? ((lj_i+4 )&15) : ((lj_i+12)&15)] <= new_lj[lj_i][6];
                        end
                    end
                end

                2'd3: begin
                    for (lj_i = 0; lj_i < CNU_COUNT; lj_i = lj_i + 1) begin
                        if (mode_reg || current_bank_sel) begin
                            Lja[0][batch ? ((lj_i+15)&15) : ((lj_i+7 )&15)] <= new_lj[lj_i][0];
                            Lja[1][batch ? ((lj_i+8 )&15) : lj_i] <= new_lj[lj_i][1];
                            Lja[2][batch ? ((lj_i+13)&15) : ((lj_i+5 )&15)] <= new_lj[lj_i][2];
                            Lja[4][batch ? ((lj_i+6 )&15) : ((lj_i+14)&15)] <= new_lj[lj_i][3];
                            Lja[5][batch ? ((lj_i+2 )&15) : ((lj_i+10)&15)] <= new_lj[lj_i][4];
                            Lja[6][batch ? ((lj_i+10)&15) : ((lj_i+2 )&15)] <= new_lj[lj_i][5];
                            Lja[7][batch ? ((lj_i+5 )&15) : ((lj_i+13)&15)] <= new_lj[lj_i][6];
                        end else begin
                            Ljb[0][batch ? ((lj_i+15)&15) : ((lj_i+7 )&15)] <= new_lj[lj_i][0];
                            Ljb[1][batch ? ((lj_i+8 )&15) : lj_i] <= new_lj[lj_i][1];
                            Ljb[2][batch ? ((lj_i+13)&15) : ((lj_i+5 )&15)] <= new_lj[lj_i][2];
                            Ljb[4][batch ? ((lj_i+6 )&15) : ((lj_i+14)&15)] <= new_lj[lj_i][3];
                            Ljb[5][batch ? ((lj_i+2 )&15) : ((lj_i+10)&15)] <= new_lj[lj_i][4];
                            Ljb[6][batch ? ((lj_i+10)&15) : ((lj_i+2 )&15)] <= new_lj[lj_i][5];
                            Ljb[7][batch ? ((lj_i+5 )&15) : ((lj_i+13)&15)] <= new_lj[lj_i][6];
                        end
                    end
                end

                default: begin
                end
            endcase
        end

        // This write is intentionally concurrent with the first LDPC batch.
        // Layer 0, batch 0 never targets column 7, position 15.
        if (late_input_write)
            Lja[7][15] <= {{2{in_data[5]}}, in_data};
    end

    integer c2v_k;
    always @(posedge clk) begin : C2V_WRITEBACK
        if (ldpc_commit) begin
            for (c2v_k = 0; c2v_k < CNU_COUNT; c2v_k = c2v_k + 1)
                c2v[layer][batch][c2v_k] <= c2v_new_lane[c2v_k];
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
