-- Copyright (C) 1999 CAL FP Bank
-- Licensed under Eiffel Forum Freeware License, version 1;
-- (see forum.txt)
--
class TEST_FFTW_OPERATOR

feature

	last_test_ok: BOOLEAN

	testit is
		do
			last_test_ok := test1
				and then test2
				and then test3
				and then test3_1
				and then test4
				and then test5

		end

	real_1d: FFTW_1D_FOURIER_OPERATOR

	in_v: ARRAY [DOUBLE] is
		once
			Result := <<1.0, 2.0, 1.0, 2.0, 1.0, 1.0, 2.0, 1.0>>
		end

feature -- 1d tests

	test1: BOOLEAN is
		local
			out_v: ARRAY [DOUBLE]
			v1: ARRAY [DOUBLE]
		do
			print ("Starting test of 1d fft... %N")
			!!out_v.make (1, 8)
			!!real_1d.make_both (8)
			real_1d.set_real (in_v)
			real_1d.set_complex (out_v)
			real_1d.forward_fft
			print ("Input after %N")
			dump (in_v)
			print ("Output after %N")			
			dump (out_v)
			print ("%N-- Done forward -- %N")
--			real_1d.set_input (out_v)
			!!v1.make (1, 8)
			real_1d.set_real (v1)
			real_1d.inverse_fft
			print ("Output of inverse %N")
			scale (v1, 1.0/ 8.0)
			dump (v1)
			print ("%N-- Done inverse --%N")
			Result := v1.is_equal (in_v)
		end


feature -- 2d tests

	real_2d: FFTW_2D_FOURIER_OPERATOR

	m_in, m_out: FFTW_MATRIX

	test2: BOOLEAN is
			-- simple forward and inverse transform test
		local
			m: FFTW_MATRIX
		do
			print ("Start 2d test...%N")
			!!real_2d.make_both (8, 8)
			!!m_in.set_spec_from_vector (in_v, 8, 8, 4, 2)
			m_in.dump
			!!m_out.set_spec_zero (8, 10)
			real_2d.set_real (m_in)
			real_2d.set_complex (m_out)
			real_2d.forward_fft
			print ("Result is: %N")
			m_out.dump

			-- try inverse
			!!m.set_spec_zero (8, 8)
--			real_2d.set_input (m_out)
			real_2d.set_real (m)
			real_2d.inverse_fft
			-- normalize result
			m.scalar_multiply (0.125*0.125)
			print ("Inverse Result is: %N")
			m.dump
			Result := m.quasi_is_equal (m_in)
		end


	test3: BOOLEAN is
			-- test linearity of the transform
		local
			m1, m2, tm1, tm2, ms, tms, tms2: FFTW_MATRIX
			v1, v2: ARRAY [DOUBLE]
			transform: FFTW_2D_FOURIER_OPERATOR
		do
			print ("Start 2d linear test %N")
			v1 := <<1.0, 2.0, 3.0, 4.0>>
			v2 := <<1.0, 1.0, 1.0, 1.0>>
			!!m1.set_spec_from_vector (v1, 4, 4, 2, 2)
			!!m2.set_spec_from_vector (v2, 4, 4, 2, 2)
			!!transform.make_both (4, 4)
			!!tm1.set_spec_zero (4, 2*((4//2)+1))
			!!tm2.set_spec_zero (4, 2*((4//2)+1))
			!!tms.set_spec_zero (4, 2*((4//2)+1))
			-- transform first matrix
			transform.set_real (m1)
			transform.set_complex (tm1)
			transform.forward_fft
			tm1.dump
			-- transform second matrix
			transform.set_real (m2)
			transform.set_complex (tm2)
			transform.forward_fft
			tm2.dump
			-- transform the sum
			ms := clone (m1)
			ms.element_add (m2)
			transform.set_real (ms)
			transform.set_complex (tms)
			transform.forward_fft
			tms.dump
			-- now see if the results agree
			tm1.element_add (tm2)
			Result := tm1.quasi_is_equal (tms)
		end


	test3_1: BOOLEAN is
			-- complete test of linearity for 2D transform
		local
			r, c: INTEGER
			transform: FFTW_2D_FOURIER_OPERATOR
			m1, m2, mc: FFTW_MATRIX
		do
			print ("Testing linearity of 2D transform...%N")
			Result := True
			!!m1.set_spec_zero (8, 8)
			!!m2.set_spec_zero (8, 2 * (8//2 + 1))
			!!transform.make_both (8, 8)
			!!mc.set_spec_zero (8, 8)
			from r := 1 until not Result or (r > 8) loop
				from c := 1 until not Result or (c > 8) loop
					m1.set_spec_zero (8, 8)
					m1.put (1.0, r, c)
					mc.set_spec_zero (8, 8)
					mc.put (1.0, r, c)
					-- forward
					transform.set_real (m1)
					transform.set_complex (m2)
					transform.forward_fft
					-- inverse
					transform.inverse_fft
					-- normalize the answer
					m1.scalar_multiply (1.0/ (8.0 * 8.0))
					Result := m1.quasi_is_equal (mc)
					if not Result then
						print ("Failed with this result: ")
						m1.dump
					end
					c := c + 1
				end
				r := r + 1
			end
		end
	


feature -- 3d tests


	test4: BOOLEAN is
			-- simple test of 3d transform
		local
			t: FFTW_3DTENSOR
			in_tensor: FFTW_3DTENSOR
			out_tensor: FFTW_3DTENSOR
			real_3d: FFTW_3D_FOURIER_OPERATOR
		do
			print ("3D test %N")
			!!in_tensor.set_spec_from_vector (in_v, 4, 4, 4, 2, 2, 2)
			!!real_3d.make_both (4, 4, 4)
			!!out_tensor.set_spec_zero (4, 4, 2 * (4//2 + 1))
			-- forward
			real_3d.set_real (in_tensor)
			real_3d.set_complex (out_tensor)
			real_3d.forward_fft
			print ("Result from forward transform %N")
			out_tensor.dump
			-- inverse
--			real_3d.set_input (out_tensor)
			!!t.set_spec_zero (4, 4, 4)
			real_3d.set_real (t)
			real_3d.inverse_fft
			print ("Result from inverse transform %N")
			t.scalar_multiply(1.0/(4.0 * 4.0 * 4.0))
			t.dump 
			Result := t.is_equal (in_tensor)
		end


	test5: BOOLEAN is
			-- test linearity of the transform
		local
			r, c, d: INTEGER
			t_in, t_out, tc: FFTW_3DTENSOR
			transform: FFTW_3D_FOURIER_OPERATOR
		do
			print ("Test linearity of 3D transform %N")
			-- transform all 4 x 4 x 4 tensors and make sure results 
			-- are as expected
			!!t_in.set_spec_zero (4, 4, 4)
			!!tc.set_spec_zero (4, 4, 4)
			!!t_out.set_spec_zero (4, 4, 2 * (4//2 + 1))
			!!transform.make_both (4, 4, 4)
			Result := True
			from r := 1 until (r > 4) or not Result loop
				from c := 1 until (c > 4) or not Result loop
					from d := 1 until (d > 4) or not Result loop
						-- prepare tensor
						t_in.set_spec_zero (4, 4, 4)
						t_in.put (1.0, r, c, d)
						tc.set_spec_zero (4, 4, 4)
						tc.put (1.0, r, c, d)
						-- forward
						transform.set_real (t_in)
						transform.set_complex (t_out)
						transform.forward_fft
						-- inverse
--						transform.set_input (t_out)
--						transform.set_output (t_in)						
						transform.inverse_fft
						t_in.scalar_multiply (1.0/(4.0*4.0*4.0))
						Result := t_in.is_equal (tc)
						if not Result then
							print ("3D failed...result: %N")
							t_in.dump
						end
						d := d + 1
					end
					c := c + 1
				end
				r := r + 1
			end
		end

	dump (v: ARRAY[DOUBLE]) is
		local
			i: INTEGER
		do
			print ("Vector size=")
			print (v.count)
			print ("%N")
			from i := 1 until i > v.count loop
				print (v @ i)
				print ("%N")
				i := i + 1
			end
			print ("------------%N")
		end

	scale (v: ARRAY [DOUBLE]; factor: DOUBLE) is
		local
			i: INTEGER
			x: DOUBLE
		do
			from i := 1 until i > v.count loop
				x := v @ i
				v.put (x * factor, i)
				i := i + 1
			end
		end

end
