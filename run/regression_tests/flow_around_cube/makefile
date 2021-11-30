.PHONY: all files build run post clean

origin_dir := ../../../
test_dir := ./build/
job_json := sgi8600.json

all: run

files:
	rsync $(origin_dir)/ $(test_dir)/ -auv \
		--exclude regression_tests/ \
		--exclude .git\* \
		--include io/input/aistwt/ --exclude io/input/*/

build: files
	make -C $(test_dir) clean
	make -C $(test_dir) -j TEST=flow_around_cube

run: build
	rsync $(job_json) $(test_dir)/run/ -v
	cd $(test_dir)/run; python3 sub.py --dir . sgi8600.json

post:
	make -C postprocess

clean:
	rm -rv $(test_dir)
	git checkout $(test_dir)
