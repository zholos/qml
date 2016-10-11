.PHONY : gen-test
gen-test:
	test/test.py -o src/test.q

bench: src_bench
