-- Copyright (C) 1999 CAL FP Bank
-- Licensed under Eiffel Forum Freeware License, version 1;
-- (see forum.txt)
--
class TEST_FFTW_MATRIX


feature -- testing matrices

	last_test_ok: BOOLEAN

	matrix: FFTW_MATRIX

	testit is
		do
			last_test_ok := test1 
				and then test2
				and then test3
				and then test4
				and then test4_1
				and then test4_2
				and then test4_3
				and then test5
				and then test6
				and then test7
				and then test8
				and then test9
				and then test10
		end


	test1: BOOLEAN is
			-- basic test of creation and setting up of the matrix
		local
			i, j: INTEGER
			k: DOUBLE
		do
			print ("Test1...%N")
			!!matrix.set_spec_zero (2, 3)
			k := 1.0
			from j := 1 until j > matrix.col_size
			loop
				from i := 1 until i > matrix.row_size
				loop
					matrix.put (k, i, j)
					k := k + 1.0
					i := i + 1
				end
				j := j + 1
			end
			matrix.dump
			-- verify that the result is correct
			Result := True
			k := 1.0
			from j := 1 until j > matrix.col_size
			loop
				from i := 1 until i > matrix.row_size
				loop
					if Result then
						Result := matrix.item (i, j) = k
					end
					k := k + 1.0
					i := i + 1
				end
				j := j + 1
			end
			print ("Test1: ")
			print (Result.out)
			print ("%N")
		end

	vect: ARRAY [DOUBLE] is
		once
			Result := <<1.0, 2.0, 3.0, 4.0>>
		end

	test2: BOOLEAN is
			-- TODO: This test should be expanded to cover bondary 
			-- conditions etc.
		local
			r, c, i: INTEGER
		do
			print ("Test2...%N")
			!!matrix.set_spec_from_vector (vect, 4, 5, 2, 2)
			matrix.dump
			Result := verify_matrix (matrix, vect, 2, 2)

			if Result then
				!!matrix.set_spec_from_vector (vect, 4, 4, 4, 1) 
				matrix.dump
				Result := verify_matrix (matrix, vect, 4, 1)
			end
			if Result then
				!!matrix.set_spec_from_vector (vect, 4, 4, 1, 4) 
				matrix.dump
				Result := verify_matrix (matrix, vect, 1, 4)
			end
			print ("Test2: ")
			print (Result)
			print ("%N")
		end

	test3: BOOLEAN is
			-- test copy and is_equal
		local
			m2: FFTW_MATRIX
		do
			print ("Test3: ... %N")
			!!matrix.set_spec_from_vector (vect, 4, 5, 1, 4)
			m2 := clone (matrix)
			m2.dump
			Result := m2.is_equal (matrix)
			print ("Test3: ")
			print (Result)
			print ("%N")
		end

	test4: BOOLEAN is
			-- test some simple operations
		local
			cp, m1, m2: FFTW_MATRIX
			r, c: INTEGER
		do
			print ("Test4: %N")
			!!m1.set_spec_from_vector (vect, 4, 4, 2, 2)
			cp := clone (m1)
			!!m2.set_spec_from_vector (vect, 4, 4, 2, 2)
			m1.real_element_multiply_by (m2)
			m1.dump
			-- check the result
			Result := True
			from c := 1 until c > m1.col_size
			loop
				from r := 1 until r > m1.row_size
				loop
					if Result then
						Result := m1.item (r, c) = (cp.item (r, c) * m2.item (r, c))
					end
					r := r + 1
				end
				c := c + 1
			end
			if Result then
				-- test multipying by a scalar
				m2.scalar_multiply (2.0)
				m2.dump
				-- verify answer
				from c := 1 until c > m2.col_size
				loop
					from r := 1 until r > m2.row_size
					loop
						if Result then
							Result := m2.item (r, c) = (cp.item (r, c) * 2.0)
						end
						r := r + 1
					end
					c := c + 1
				end
				
			end

			print ("Test 4: ")
			print (Result)
			print ("%N")
		end

	test4_1: BOOLEAN is
			-- test conversions
		local
			m: FFTW_MATRIX
			v: ARRAY [DOUBLE]
		do
			print ("Testing conversions to vector...%N")
			!!m.set_spec_from_vector (vect, 3, 3, 2, 2)
			v := m.to_vector (2, 2)
			Result := v.is_equal(vect)
		end

	test4_2: BOOLEAN is
			-- test conjugation
		local
			m1, m2: FFTW_MATRIX
		do
			print ("Testing congugation %N")
			!!m1.set_spec_from_vector (vect, 3, 6, 2, 2)
			m2 := clone (m1)
			print ("Before %N")
			m2.dump
			m1.conjugate
			print ("After %N")
			m1.dump
			Result := not m1.is_equal (m2)
			if Result then
				-- congugate twice, should be equal
				m1.conjugate
				Result := m1.is_equal (m2)
			end
		end

	test4_3: BOOLEAN is
			-- test complex multiply
		local
			v: ARRAY [DOUBLE]
			m1, m2: FFTW_MATRIX
		do
			print ("Test complex multiply %N")
			v := <<0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6>>
			!!m1.set_spec_from_vector (v, 3, 4, 3, 4)
			m2 := clone (m1)
			m1.complex_element_multiply_by (m2)
			print ("Result: %N")
			m1.dump
			-- check result
			v := <<-1, 0, -4, 0, -9, 0, -16, 0, -25, 0, -36, 0>>
			!!m2.set_spec_from_vector (v, 3, 4, 3, 4)
			Result := m1.is_equal (m2)
		end


feature -- testing tensors

	tensor: FFTW_3DTENSOR

	t_vect: ARRAY [DOUBLE] is
		once
			Result := <<1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 , 8.0>>
		end

	test5: BOOLEAN is
			-- basic creation test
		local
			v: ARRAY [DOUBLE]
		do
			print ("Testing 3d tensor%N")
			!!tensor.set_spec_from_vector (t_vect, 3, 3, 3, 2, 2, 2)
			tensor.dump
			-- verify the result
			Result := verify_tensor (tensor, t_vect, 2, 2, 2)
		end


	test6: BOOLEAN is
			-- copy and is_equal
		local
			t1, t2: FFTW_3DTENSOR
		do
			print ("Testing copying and comparison of 3d tensors %N")
			!!t1.set_spec_from_vector (t_vect, 4, 4, 4, 4, 2, 1)
			t1.dump
			t2 := clone (t1)
			Result := t2.is_equal (t1)
		end


	test7: BOOLEAN is
			-- test multiplication
		local
			t1, t2, t3: FFTW_3DTENSOR
			v: ARRAY [DOUBLE]
			r, c, d: INTEGER
		do
			print ("Testing multiplication by a scalar...%N")
			!!t1.set_spec_from_vector (t_vect, 3, 3, 3, 2, 2, 2)
			t2 := clone (t1)
			t1.scalar_multiply (2.0)
			t1.dump
			t1.scalar_multiply (0.5)
			Result := t1.is_equal (t2)
			if not Result then
				print ("scalar_multiply failed...%N")
			else
				print ("Testing elemnt-wise multiplication %N")
				-- Try multiplying one by another
				t3 := clone (t1)
				t1.real_element_multiply_by (t2)
				t1.dump
				-- verify the result
				from r := 1 until r > t1.row_size loop
					from c := 1 until c > t1.col_size loop
						from d := 1 until d > t1.depth_size loop
							if Result then
								Result := t1.item (r, c, d) = t2.item (r, c, d) * t3.item (r, c, d)
							end
							d := d + 1
						end
						c := c + 1
					end
					r := r + 1
				end
					
			end
		end

	test8: BOOLEAN is
		local
			v: ARRAY [DOUBLE]
			t1: FFTW_3DTENSOR
		do
			print ("Test conversions...%N")
			!!t1.set_spec_from_vector (t_vect, 3, 3, 3, 2, 2, 2)
			v := t1.to_vector (2, 2, 2)
			Result := v.is_equal (t_vect)
		end


	test9: BOOLEAN is
			-- test conjugation
		local
			t1, t2: FFTW_3DTENSOR
		do
			print ("Test tensor conjugation...%N")
			!!t1.set_spec_from_vector (t_vect, 4, 4, 4, 2, 2, 2)
			t2 := clone (t1)
			t1.conjugate
			t1.dump
			Result := not t1.is_equal (t2)
			if Result then
				t1.conjugate
				Result := t1.is_equal (t2)
			end
		end
			

	test10: BOOLEAN is
			-- test complex multiply
		local
			t1, t2, t3: FFTW_3DTENSOR
			v: ARRAY [DOUBLE]
		do
			print ("Test complex multiply %N")
			v := <<0, 1, 0, 1, 0, 1, 0, 1>>
			!!t1.set_spec_from_vector (v, 2, 2, 2, 2, 2, 2)
			t2 := clone (t1)
			t1.complex_element_multiply_by (t2)
			t1.dump
			v := <<-1, 0, -1, 0, -1, 0, -1, 0>>
			!!t3.set_spec_from_vector (v, 2, 2, 2, 2, 2, 2)
			Result := t1.is_equal (t3)
		end

feature -- helper routines

	verify_matrix (m: FFTW_MATRIX; v: ARRAY [DOUBLE]; lr, lc: INTEGER): BOOLEAN is
		local
			i, c, r: INTEGER
		do			
			-- verify that the result is correct
			Result := True
			i := 1
			from r := 1 until r > m.row_size
			loop
				from c := 1 until c > m.col_size
				loop
					if Result then
						if (c <= lc) and (r <= lr) then
							Result := (v @ i) = m.item (r, c)
							i := i + 1
						else
							Result  := m.item (r, c) = 0.0
						end
					end
					c := c + 1
				end
				r := r + 1
			end
		end


	verify_tensor (t: FFTW_3DTENSOR; v: ARRAY [DOUBLE]; lr, lc, ld: INTEGER): BOOLEAN is
		local
			i, r, c, d: INTEGER
		do
			Result := True
			i := 1
			from r := 1 until r > t.row_size loop
				from c := 1 until c > t.col_size loop
					from d := 1 until d > t.depth_size loop
						if Result then
							if (r <= lr) and (c <= lc) and (d <= ld) then
								Result :=  v.item (i) = t.item (r, c, d)
								i := i + 1
							else
								Result :=  t.item (r, c, d) = 0.0
							end
						end
						d := d + 1
					end
					c := c + 1
				end
				r := r + 1
			end

		end

end
