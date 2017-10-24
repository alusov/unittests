all dep clean indent tests::
	cd gtest/lib && make $@ && cd .. \\
	cd testfuncs && make $@ && cd .. \\
	cd testintervals && make $@ && cd .. \\
	cd testglobopt && make $@ && cd .. \\
	cd testder && make $@ && cd .. \\

