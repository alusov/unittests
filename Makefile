all dep clean indent tests::
	cd testfuncs && make $@ && cd .. \\
	cd testintervals && make $@ && cd ..
