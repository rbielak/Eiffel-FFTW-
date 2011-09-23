-- Copyright (C) 1999 CAL FP Bank
-- Licensed under Eiffel Forum Freeware License, version 1;
-- (see forum.txt)
--
class TEST_ALL

creation

	make

feature

	test1: TEST_FFTW_MATRIX

	test2: TEST_FFTW_OPERATOR

	make is
		do
			!!test1
			!!test2
			test1.testit
			if not test1.last_test_ok then
				print ("*** Test1 failed %N")
			else
				test2.testit
				if not test2.last_test_ok then
					print ("*** Test2 failed %N")
				end
			end
			print ("*** Done... %N")
		end

end
