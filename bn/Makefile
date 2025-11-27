BUILD_DIR = build
INC_DIR = inc
SRC_DIR = src

all:
	gcc \
		-Wall -Wextra \
		-O3 -funroll-loops \
		$(SRC_DIR)/common.c \
		$(SRC_DIR)/wlc_bench.c  \
		$(SRC_DIR)/wlc_bn.c \
		$(SRC_DIR)/wlc_rand_stdc.c \
		$(SRC_DIR)/wlc_test.c \
		$(SRC_DIR)/wlc_util.c \
        $(SRC_DIR)/wlc_bn_print.c  \
		$(SRC_DIR)/main.c \
		-o $(BUILD_DIR)/main.exe


debug:
	gcc \
                -Wall -Wextra \
                -g \
                $(SRC_DIR)/common.c \
                $(SRC_DIR)/wlc_bench.c  \
                $(SRC_DIR)/wlc_bn.c \
                $(SRC_DIR)/wlc_rand_stdc.c \
                $(SRC_DIR)/wlc_test.c \
                $(SRC_DIR)/wlc_util.c \
                $(SRC_DIR)/wlc_bn_print.c  \
                $(SRC_DIR)/main.c \
                -o $(BUILD_DIR)/main.exe


test:
	        gcc \
                -Wall -Wextra \
                -O3 -funroll-loops \
                $(SRC_DIR)/common.c \
                $(SRC_DIR)/wlc_bench.c  \
                $(SRC_DIR)/wlc_bn.c \
                $(SRC_DIR)/wlc_rand_stdc.c \
                $(SRC_DIR)/wlc_test.c \
                $(SRC_DIR)/wlc_util.c \
                $(SRC_DIR)/wlc_bn_print.c  \
                $(SRC_DIR)/test.c \
                -o $(BUILD_DIR)/test.exe


bench:
	        gcc \
                -Wall -Wextra \
                -O3 -funroll-loops \
                $(SRC_DIR)/common.c \
                $(SRC_DIR)/wlc_bench.c  \
                $(SRC_DIR)/wlc_bn.c \
                $(SRC_DIR)/wlc_rand_stdc.c \
                $(SRC_DIR)/wlc_test.c \
                $(SRC_DIR)/wlc_util.c \
                $(SRC_DIR)/wlc_bn_print.c  \
                $(SRC_DIR)/bench.c \
                -o $(BUILD_DIR)/bench.exe


clean:
	rm -f $(BUILD_DIR)/main.exe $(BUILD_DIR)/test.exe $(BUILD_DIR)/bench.exe
