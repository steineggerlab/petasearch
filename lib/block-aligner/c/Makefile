CC=clang
CFLAGS=-L../target/release -lblock_aligner -I.

.PHONY: all block_aligner debug

all: block_aligner debug

block_aligner:
	BLOCK_ALIGNER_C= cargo build --release --features simd_avx2

debug:
	RUSTFLAGS=-g cargo build --features simd_avx2
