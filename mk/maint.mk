.PHONY : gen-test
gen-test:
	test/test.py -o src/test.q

bench: src_bench

git-release:
	git checkout -b release
	sed -i '/^VERSION :=/s/-.*//' mk/src.mk
	git add mk/src.mk
	git commit -m 'fixed version number'
	git checkout master
	git merge -X theirs --no-edit release
	git branch -d release
	VERSION=`sed -n '/^VERSION := */s///p' mk/src.mk`; \
	    git tag -a "v$$VERSION" -m "version $$VERSION"
