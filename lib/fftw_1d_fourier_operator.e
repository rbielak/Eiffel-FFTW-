-- Copyright (C) 1999 CAL FP Bank
-- Licensed under Eiffel Forum Freeware License, version 1;
-- (see forum.txt)
--
indexing

	description: "1 dimensional real fourier operator, using FFTW. %
                 %The complete reference for the FFTW code can be %
                 %found at http://www.fftw.org"
	 author: "Richie Bielak"
	 date: "December 10th, 1999"

class FFTW_1D_FOURIER_OPERATOR

inherit

	FOURIER_OPERATOR

	FFTW_EXTERNALS

	MEMORY
		redefine
			dispose
		end

creation

	make_both,
	make_forward_plan_only,
	make_inverse_plan_only
	

feature -- creation


	make_both (lsize: INTEGER) is
			-- create plans for both forward and inverse transformations
		require
			valid_size: lsize > 0
		do
			make_forward_plan_only (lsize)
			make_inverse_plan_only (lsize)
		ensure
			size_set: size = lsize
			iplan_valid: inverse_plan /= default_pointer
			fplan_valid: forward_plan /= default_pointer
		end

	make_forward_plan_only  (lsize: INTEGER) is
			-- create plans for forward transformation only
		require
			valid_size: lsize > 0
		do
			size := lsize
			-- free previous plan if exists
			if forward_plan /= default_pointer then
				rfftw_destroy_plan (forward_plan)
			end
			forward_plan := rfftw_create_plan (size, -1, 0)
		ensure
			size_set: size = lsize
			plan_valid: forward_plan /= default_pointer
		end

	make_inverse_plan_only  (lsize: INTEGER) is
			-- create plans for forward transformation only
		require
			valid_size: lsize > 0
		do
			size := lsize
			if inverse_plan /= default_pointer then
				rfftw_destroy_plan (inverse_plan)
			end
			inverse_plan := rfftw_create_plan (size, 1, 0) 
		ensure
			size_set: size = lsize
			plan_valid: inverse_plan /= default_pointer
		end

feature -- input and output

	size: INTEGER
			-- size of input vector. The transform will be fastest 
			-- when this number is a power of 2.

	real_vector: ARRAY [DOUBLE]
			-- this is the input for the forward transform and output 
			-- of the inverse transform

	set_real (v: like real_vector) is
		require
			valid_vector: (v /= Void) and then ((v.count = size) and (v.lower = 1))
		do
			real_vector := v
		ensure
			set: real_vector = v
		end

	complex_vector: like real_vector
			-- this is the output of forward fft and input to inverse fft.
			-- This vector is  halfcomplex: a halfcomplex vector contains the real 
			-- and imaginary parts. For each number the real part is 
			-- in v[i] and imaginary in v[count - i], except for v[1] 
			-- and v[count/2]. 
	
	
	set_complex (v: like complex_vector) is
		require
			valid_vector: (v /= Void) and then ((v.count = size) and (v.lower = 1))
		do
			complex_vector := v
		ensure
			set: complex_vector = v
		end

feature -- operations

	forward_fft is
		do
			rfftw_one (forward_plan, $(real_vector.to_c), $(complex_vector.to_c))
			last_direction := 1
		end
		
	inverse_fft is
		do
			rfftw_one (inverse_plan, $(complex_vector.to_c), $(real_vector.to_c))
			last_direction := -1
		end

feature -- queries

	
	ready_for_forward: BOOLEAN is
		do
			Result := (forward_plan /= default_pointer) 
				and (real_vector /= Void)
				and (complex_vector /= Void)
		end

	ready_for_inverse: BOOLEAN is
		do
			Result := (inverse_plan /= default_pointer)
				and (real_vector /= Void)
				and (complex_vector /= Void)
		end

feature {NONE} -- implementation

	forward_plan: POINTER

	inverse_plan: POINTER

	dispose is
		do
			if forward_plan /= default_pointer then
				rfftw_destroy_plan (forward_plan)
			end
			if inverse_plan /= default_pointer then
				rfftw_destroy_plan (inverse_plan)
			end
		end


invariant

	plan_defined: (forward_plan /= default_pointer) or (inverse_plan /= default_pointer)
	consistent_dimension:  size > 0

end
