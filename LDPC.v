module cnu7_desc (
    input      [41:0] q_flat,
    output reg [19:0] desc_out,
    output reg [41:0] r_flat
);

    reg signed [5:0] q_value [0:6];
    reg        [4:0] q_mag   [0:6];
    reg        [6:0] q_neg;
    reg        [6:0] r_neg;
    reg              sign_parity;
    reg        [12:0] leaf0;
    reg        [12:0] leaf1;
    reg        [12:0] leaf2;
    reg        [12:0] leaf3;
    reg        [12:0] leaf4;
    reg        [12:0] leaf5;
    reg        [12:0] leaf6;
    reg        [12:0] pair01;
    reg        [12:0] pair23;
    reg        [12:0] pair45;
    reg        [12:0] group03;
    reg        [12:0] group46;
    reg        [12:0] root_pair;
    reg        [4:0] norm_a;
    reg        [4:0] norm_b;
    reg        [4:0] out_mag;
    integer edge_no;

    function [12:0] merge_pair;
        input [12:0] left_value;
        input [12:0] right_value;
        reg [4:0] left_min;
        reg [4:0] left_second;
        reg [2:0] left_index;
        reg [4:0] right_min;
        reg [4:0] right_second;
        reg [2:0] right_index;
        reg [4:0] next_second;
        begin
            left_min    = left_value[12:8];
            left_index  = left_value[7:5];
            left_second = left_value[4:0];
            right_min    = right_value[12:8];
            right_index  = right_value[7:5];
            right_second = right_value[4:0];

            if ((left_min < right_min) ||
                ((left_min == right_min) && (left_index <= right_index))) begin
                if (left_second < right_min)
                    next_second = left_second;
                else
                    next_second = right_min;
                merge_pair = {left_min, left_index, next_second};
            end
            else begin
                if (right_second < left_min)
                    next_second = right_second;
                else
                    next_second = left_min;
                merge_pair = {right_min, right_index, next_second};
            end
        end
    endfunction

    function [4:0] norm5;
        input [4:0] mag;
        reg [1:0] rem_round;
        begin
            case (mag[1:0])
                2'd0: rem_round = 2'd0;
                2'd1: rem_round = 2'd1;
                2'd2: rem_round = 2'd2;
                default: rem_round = 2'd2;
            endcase
            norm5 = {2'b00, mag[4:2]} +
                    {1'b0, mag[4:2], 1'b0} +
                    {3'b000, rem_round};
        end
    endfunction

    always @* begin
        q_neg = 7'd0;
        r_neg = 7'd0;
        r_flat = 42'd0;
        desc_out = 20'd0;

        for (edge_no = 0; edge_no < 7; edge_no = edge_no + 1) begin
            q_value[edge_no] = q_flat[edge_no*6 +: 6];
            q_neg[edge_no] = q_value[edge_no][5];
            if (q_value[edge_no][5])
                q_mag[edge_no] = (~q_value[edge_no][4:0]) + 5'd1;
            else
                q_mag[edge_no] = q_value[edge_no][4:0];
        end

        leaf0 = {q_mag[0], 3'd0, 5'd31};
        leaf1 = {q_mag[1], 3'd1, 5'd31};
        leaf2 = {q_mag[2], 3'd2, 5'd31};
        leaf3 = {q_mag[3], 3'd3, 5'd31};
        leaf4 = {q_mag[4], 3'd4, 5'd31};
        leaf5 = {q_mag[5], 3'd5, 5'd31};
        leaf6 = {q_mag[6], 3'd6, 5'd31};

        pair01 = merge_pair(leaf0, leaf1);
        pair23 = merge_pair(leaf2, leaf3);
        pair45 = merge_pair(leaf4, leaf5);
        group03 = merge_pair(pair01, pair23);
        group46 = merge_pair(pair45, leaf6);
        root_pair = merge_pair(group03, group46);

        norm_a = norm5(root_pair[12:8]);
        norm_b = norm5(root_pair[4:0]);
        sign_parity = ^q_neg;

        for (edge_no = 0; edge_no < 7; edge_no = edge_no + 1) begin
            r_neg[edge_no] = sign_parity ^ q_neg[edge_no];
            if (edge_no[2:0] == root_pair[7:5])
                out_mag = norm_b;
            else
                out_mag = norm_a;

            if (out_mag == 5'd0)
                r_flat[edge_no*6 +: 6] = 6'd0;
            else if (r_neg[edge_no])
                r_flat[edge_no*6 +: 6] = (~{1'b0, out_mag}) + 6'd1;
            else
                r_flat[edge_no*6 +: 6] = {1'b0, out_mag};
        end

        desc_out = {norm_a, norm_b, root_pair[7:5], r_neg};
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

    localparam [2:0] ST_IDLE   = 3'd0;
    localparam [2:0] ST_LOAD   = 3'd1;
    localparam [2:0] ST_RUN    = 3'd2;
    localparam [2:0] ST_CHECK  = 3'd3;
    localparam [2:0] ST_OUTPUT = 3'd4;

    localparam integer CNU_COUNT = 8;
    localparam integer EDGE_COUNT = 56;

    reg [2:0] state;
    reg       mode_reg;
    reg [6:0] input_idx;
    reg [6:0] output_idx;
    reg [2:0] iter_idx;
    reg [1:0] layer_cnt;
    reg [1:0] batch_cnt;
    reg       bank_sel;
    reg       warn_reg;

    reg signed [7:0] p_mem0 [0:127];
    reg signed [7:0] p_mem1 [0:127];
    reg        [19:0] desc_mem [0:63];

    reg  [CNU_COUNT*42-1:0] q_bus;
    wire [CNU_COUNT*42-1:0] r_bus;
    wire [CNU_COUNT*20-1:0] desc_bus;
    reg  [CNU_COUNT*6-1:0]  cn_idx_bus;
    reg  [EDGE_COUNT*7-1:0] vn_idx_bus;
    reg  [EDGE_COUNT*8-1:0] q_base_bus;
    reg  [EDGE_COUNT*8-1:0] acc_base_bus;
    reg  [EDGE_COUNT*8-1:0] p_new_bus;

    reg [63:0] syn_bits;
    reg [6:0] syn_inputs;
    reg       syn_nonzero;

    wire [1:0] phase_layer;
    wire [1:0] phase_batch;
    wire       zero_old;
    wire       at_limit;
    wire       advance_iter;
    wire       batch_commit;
    wire       phase_active;

    assign phase_layer = (state == ST_CHECK) ? 2'd0 : layer_cnt;
    assign phase_batch = (state == ST_CHECK) ? 2'd0 : batch_cnt;
    assign zero_old = (state == ST_RUN) && (iter_idx == 3'd0);
    assign at_limit = (iter_idx == 3'd7);
    assign advance_iter = syn_nonzero && !at_limit;
    assign batch_commit = (state == ST_RUN) ||
                          ((state == ST_CHECK) && advance_iter);
    assign phase_active = (state == ST_RUN) || (state == ST_CHECK);

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

    function [2:0] edge_column;
        input [1:0] layer_value;
        input [2:0] edge_value;
        begin
            if (edge_value < {1'b0, layer_value})
                edge_column = edge_value;
            else
                edge_column = edge_value + 3'd1;
        end
    endfunction

    function [3:0] shift_value;
        input [1:0] layer_value;
        input [2:0] column_value;
        begin
            shift_value = 4'd0;
            case (layer_value)
                2'd0: begin
                    case (column_value)
                        3'd1: shift_value = 4'd14;
                        3'd2: shift_value = 4'd10;
                        3'd3: shift_value = 4'd2;
                        3'd4: shift_value = 4'd13;
                        3'd5: shift_value = 4'd12;
                        3'd6: shift_value = 4'd9;
                        3'd7: shift_value = 4'd3;
                        default: shift_value = 4'd0;
                    endcase
                end
                2'd1: begin
                    case (column_value)
                        3'd0: shift_value = 4'd5;
                        3'd2: shift_value = 4'd14;
                        3'd3: shift_value = 4'd10;
                        3'd4: shift_value = 4'd2;
                        3'd5: shift_value = 4'd13;
                        3'd6: shift_value = 4'd12;
                        3'd7: shift_value = 4'd9;
                        default: shift_value = 4'd0;
                    endcase
                end
                2'd2: begin
                    case (column_value)
                        3'd0: shift_value = 4'd0;
                        3'd1: shift_value = 4'd5;
                        3'd3: shift_value = 4'd14;
                        3'd4: shift_value = 4'd10;
                        3'd5: shift_value = 4'd2;
                        3'd6: shift_value = 4'd13;
                        3'd7: shift_value = 4'd12;
                        default: shift_value = 4'd0;
                    endcase
                end
                2'd3: begin
                    case (column_value)
                        3'd0: shift_value = 4'd7;
                        3'd1: shift_value = 4'd0;
                        3'd2: shift_value = 4'd5;
                        3'd4: shift_value = 4'd14;
                        3'd5: shift_value = 4'd10;
                        3'd6: shift_value = 4'd2;
                        3'd7: shift_value = 4'd13;
                        default: shift_value = 4'd0;
                    endcase
                end
                default: shift_value = 4'd0;
            endcase
        end
    endfunction

    genvar gen_lane;
    generate
        for (gen_lane = 0; gen_lane < CNU_COUNT; gen_lane = gen_lane + 1) begin : GEN_CNU
            cnu7_desc u_cnu7_desc (
                .q_flat  (q_bus[gen_lane*42 +: 42]),
                .desc_out(desc_bus[gen_lane*20 +: 20]),
                .r_flat  (r_bus[gen_lane*42 +: 42])
            );
        end
    endgenerate

    integer lane_no;
    integer route_edge;
    reg [4:0] route_row;
    reg [5:0] route_cn;
    reg [2:0] route_col;
    reg [3:0] route_shift;
    reg [3:0] route_x;
    reg [6:0] route_vn;
    reg [19:0] route_desc;
    reg [4:0] old_mag;
    reg signed [5:0] old_r;
    reg signed [7:0] old_r_ext;
    reg signed [7:0] source_value;
    reg signed [7:0] accum_value;
    reg signed [7:0] q_base_value;
    reg signed [7:0] acc_base_value;

    always @* begin
        q_bus = {(CNU_COUNT*42){1'b0}};
        cn_idx_bus = {(CNU_COUNT*6){1'b0}};
        vn_idx_bus = {(EDGE_COUNT*7){1'b0}};
        q_base_bus = {(EDGE_COUNT*8){1'b0}};
        acc_base_bus = {(EDGE_COUNT*8){1'b0}};

        route_row = 5'd0;
        route_cn = 6'd0;
        route_col = 3'd0;
        route_shift = 4'd0;
        route_x = 4'd0;
        route_vn = 7'd0;
        route_desc = 20'd0;
        old_mag = 5'd0;
        old_r = 6'sd0;
        old_r_ext = 8'sd0;
        source_value = 8'sd0;
        accum_value = 8'sd0;
        q_base_value = 8'sd0;
        acc_base_value = 8'sd0;

        if (phase_active) begin
            for (lane_no = 0; lane_no < CNU_COUNT; lane_no = lane_no + 1) begin
                route_row = {phase_batch, 3'b000} + lane_no[4:0];
                route_cn = {phase_layer, 4'b0000} + route_row;
                cn_idx_bus[lane_no*6 +: 6] = route_cn;

                if (zero_old)
                    route_desc = 20'd0;
                else
                    route_desc = desc_mem[route_cn];

                for (route_edge = 0; route_edge < 7; route_edge = route_edge + 1) begin
                    route_col = edge_column(phase_layer, route_edge[2:0]);
                    route_shift = shift_value(phase_layer, route_col);
                    route_x = route_row[3:0] + route_shift;
                    route_vn = {route_col, 4'b0000} + {3'b000, route_x};
                    vn_idx_bus[(lane_no*7+route_edge)*7 +: 7] = route_vn;

                    if (route_edge[2:0] == route_desc[9:7])
                        old_mag = route_desc[14:10];
                    else
                        old_mag = route_desc[19:15];

                    if (zero_old || (old_mag == 5'd0))
                        old_r = 6'sd0;
                    else if (route_desc[route_edge])
                        old_r = (~{1'b0, old_mag}) + 6'd1;
                    else
                        old_r = {1'b0, old_mag};

                    old_r_ext = {{2{old_r[5]}}, old_r};

                    if (mode_reg)
                        source_value = p_mem0[route_vn];
                    else if (bank_sel)
                        source_value = p_mem1[route_vn];
                    else
                        source_value = p_mem0[route_vn];

                    if (mode_reg)
                        accum_value = source_value;
                    else if (((route_col == 3'd0) && (phase_layer == 2'd1)) ||
                             ((route_col != 3'd0) && (phase_layer == 2'd0)))
                        accum_value = source_value;
                    else if (bank_sel)
                        accum_value = p_mem0[route_vn];
                    else
                        accum_value = p_mem1[route_vn];

                    q_base_value = source_value - old_r_ext;
                    acc_base_value = accum_value - old_r_ext;
                    q_bus[(lane_no*42+route_edge*6) +: 6] = clip6(q_base_value);
                    q_base_bus[(lane_no*56+route_edge*8) +: 8] = q_base_value;
                    acc_base_bus[(lane_no*56+route_edge*8) +: 8] = acc_base_value;
                end
            end
        end
    end

    integer post_lane;
    integer post_edge;
    reg signed [5:0] post_r;
    reg signed [7:0] post_r_ext;
    reg signed [7:0] post_q_base;
    reg signed [7:0] post_acc_base;
    reg signed [7:0] post_update_base;
    reg signed [7:0] post_value;

    always @* begin
        p_new_bus = {(EDGE_COUNT*8){1'b0}};
        post_r = 6'sd0;
        post_r_ext = 8'sd0;
        post_q_base = 8'sd0;
        post_acc_base = 8'sd0;
        post_update_base = 8'sd0;
        post_value = 8'sd0;

        for (post_lane = 0; post_lane < CNU_COUNT; post_lane = post_lane + 1) begin
            for (post_edge = 0; post_edge < 7; post_edge = post_edge + 1) begin
                post_r = r_bus[(post_lane*42+post_edge*6) +: 6];
                post_r_ext = {{2{post_r[5]}}, post_r};
                post_q_base = q_base_bus[(post_lane*56+post_edge*8) +: 8];
                post_acc_base = acc_base_bus[(post_lane*56+post_edge*8) +: 8];
                if (mode_reg)
                    post_update_base = post_q_base;
                else
                    post_update_base = post_acc_base;
                post_value = post_update_base + post_r_ext;
                p_new_bus[(post_lane*56+post_edge*8) +: 8] = post_value;
            end
        end
    end

    integer syn_layer;
    integer syn_row;
    integer syn_edge;
    reg [2:0] syn_col;
    reg [3:0] syn_shift;
    reg [3:0] syn_x;
    reg [6:0] syn_vn;

    always @* begin
        syn_bits = 64'd0;
        syn_inputs = 7'd0;
        syn_col = 3'd0;
        syn_shift = 4'd0;
        syn_x = 4'd0;
        syn_vn = 7'd0;

        for (syn_layer = 0; syn_layer < 4; syn_layer = syn_layer + 1) begin
            for (syn_row = 0; syn_row < 16; syn_row = syn_row + 1) begin
                for (syn_edge = 0; syn_edge < 7; syn_edge = syn_edge + 1) begin
                    syn_col = edge_column(syn_layer[1:0], syn_edge[2:0]);
                    syn_shift = shift_value(syn_layer[1:0], syn_col);
                    syn_x = syn_row[3:0] + syn_shift;
                    syn_vn = {syn_col, 4'b0000} + {3'b000, syn_x};
                    if (mode_reg)
                        syn_inputs[syn_edge] = p_mem0[syn_vn][7];
                    else if (bank_sel)
                        syn_inputs[syn_edge] = p_mem1[syn_vn][7];
                    else
                        syn_inputs[syn_edge] = p_mem0[syn_vn][7];
                end
                syn_bits[syn_layer*16+syn_row] =
                    (syn_inputs[0] ^ syn_inputs[1]) ^
                    (syn_inputs[2] ^ syn_inputs[3]) ^
                    (syn_inputs[4] ^ syn_inputs[5]) ^ syn_inputs[6];
            end
        end
        syn_nonzero = |syn_bits;
    end

    integer write_lane;
    integer write_edge;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state      <= ST_IDLE;
            mode_reg   <= 1'b0;
            input_idx  <= 7'd0;
            output_idx <= 7'd0;
            iter_idx   <= 3'd0;
            layer_cnt  <= 2'd0;
            batch_cnt  <= 2'd0;
            bank_sel   <= 1'b0;
            warn_reg   <= 1'b0;
            out_valid  <= 1'b0;
            out_data   <= 8'd0;
            out_warn   <= 1'b0;
        end
        else begin
            if (batch_commit) begin
                for (write_lane = 0; write_lane < CNU_COUNT; write_lane = write_lane + 1) begin
                    desc_mem[cn_idx_bus[write_lane*6 +: 6]] <=
                        desc_bus[write_lane*20 +: 20];
                    for (write_edge = 0; write_edge < 7; write_edge = write_edge + 1) begin
                        if (mode_reg)
                            p_mem0[vn_idx_bus[(write_lane*7+write_edge)*7 +: 7]] <=
                                p_new_bus[(write_lane*56+write_edge*8) +: 8];
                        else if (bank_sel)
                            p_mem0[vn_idx_bus[(write_lane*7+write_edge)*7 +: 7]] <=
                                p_new_bus[(write_lane*56+write_edge*8) +: 8];
                        else
                            p_mem1[vn_idx_bus[(write_lane*7+write_edge)*7 +: 7]] <=
                                p_new_bus[(write_lane*56+write_edge*8) +: 8];
                    end
                end
            end

            case (state)
                ST_IDLE: begin
                    out_valid <= 1'b0;
                    out_data <= 8'd0;
                    out_warn <= 1'b0;
                    warn_reg <= 1'b0;
                    if (in_mode_valid) begin
                        mode_reg <= in_mode;
                        input_idx <= 7'd0;
                        output_idx <= 7'd0;
                        iter_idx <= 3'd0;
                        layer_cnt <= 2'd0;
                        batch_cnt <= 2'd0;
                        bank_sel <= 1'b0;
                        state <= ST_LOAD;
                    end
                end

                ST_LOAD: begin
                    out_valid <= 1'b0;
                    out_data <= 8'd0;
                    out_warn <= 1'b0;
                    if (in_data_valid) begin
                        p_mem0[input_idx] <= {{2{in_data[5]}}, in_data};
                        if (input_idx == 7'd127) begin
                            input_idx <= 7'd0;
                            iter_idx <= 3'd0;
                            layer_cnt <= 2'd0;
                            batch_cnt <= 2'd0;
                            bank_sel <= 1'b0;
                            state <= ST_RUN;
                        end
                        else begin
                            input_idx <= input_idx + 7'd1;
                        end
                    end
                end

                ST_RUN: begin
                    out_valid <= 1'b0;
                    out_data <= 8'd0;
                    out_warn <= 1'b0;
                    if (batch_cnt == 2'd1) begin
                        batch_cnt <= 2'd0;
                        if (layer_cnt == 2'd3) begin
                            layer_cnt <= 2'd0;
                            if (!mode_reg)
                                bank_sel <= ~bank_sel;
                            state <= ST_CHECK;
                        end
                        else begin
                            layer_cnt <= layer_cnt + 2'd1;
                        end
                    end
                    else begin
                        batch_cnt <= batch_cnt + 2'd1;
                    end
                end

                ST_CHECK: begin
                    if (advance_iter) begin
                        out_valid <= 1'b0;
                        out_data <= 8'd0;
                        out_warn <= 1'b0;
                        iter_idx <= iter_idx + 3'd1;
                        layer_cnt <= 2'd0;
                        batch_cnt <= 2'd1;
                        state <= ST_RUN;
                    end
                    else begin
                        warn_reg <= syn_nonzero && at_limit;
                        out_valid <= 1'b1;
                        if (mode_reg)
                            out_data <= p_mem0[0];
                        else if (bank_sel)
                            out_data <= p_mem1[0];
                        else
                            out_data <= p_mem0[0];
                        out_warn <= syn_nonzero && at_limit;
                        output_idx <= 7'd1;
                        state <= ST_OUTPUT;
                    end
                end

                ST_OUTPUT: begin
                    out_valid <= 1'b1;
                    out_warn <= warn_reg;
                    if (mode_reg)
                        out_data <= p_mem0[output_idx];
                    else if (bank_sel)
                        out_data <= p_mem1[output_idx];
                    else
                        out_data <= p_mem0[output_idx];

                    if (output_idx == 7'd127) begin
                        output_idx <= 7'd0;
                        state <= ST_IDLE;
                    end
                    else begin
                        output_idx <= output_idx + 7'd1;
                    end
                end

                default: begin
                    state <= ST_IDLE;
                    mode_reg <= 1'b0;
                    input_idx <= 7'd0;
                    output_idx <= 7'd0;
                    iter_idx <= 3'd0;
                    layer_cnt <= 2'd0;
                    batch_cnt <= 2'd0;
                    bank_sel <= 1'b0;
                    warn_reg <= 1'b0;
                    out_valid <= 1'b0;
                    out_data <= 8'd0;
                    out_warn <= 1'b0;
                end
            endcase
        end
    end

endmodule
