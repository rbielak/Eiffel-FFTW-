-- Copyright (C) 1999 CAL FP Bank
-- Licensed under Eiffel Forum Freeware License, version 1;
-- (see forum.txt)
--
indexing

	description: "Eiffel representation of a 3D tensor. These things %
                 %are inputs into FFT"

class FFTW_3DTENSOR

inherit

	EXCEPTIONS
		redefine
			copy, is_equal
		end

creation

	set_spec_from_vector, 
	set_spec_zero,
	make

feature -- creation

	set_spec_from_vector (v: ARRAY [DOUBLE]; big_row, big_col, big_depth, little_row, little_col, little_depth: INTEGER) is
		require
			valid_vector: v /= Void
			consistent_row: (big_row >= little_row) and (little_row > 0)
			consistent_column: (big_col >= little_col) and (little_col > 0)
			consistent_depth: (big_depth >= little_depth) and (little_depth > 0)
			consistent_with_vector: (little_col * little_row * little_depth) = v.count
		local
			i, r, c, d: INTEGER
			to_area: TO_SPECIAL [DOUBLE]
		do
			make (big_row, big_col, big_depth)
			i := 1
			from r := 1 until r > row_size loop
				from c := 1 until c > col_size loop
					from d := 1 until d > depth_size loop
						if (r <= little_row) and (c <= little_col) and (d <= little_depth) then
							-- put (v @ i, r, c, d)
							area.put (v @ i, (d - 1) + depth_size * ((c - 1) + (r - 1) * col_size))
							i := i + 1
						else
							-- put (0.0, r, c, d)
							area.put (0.0, (d - 1) + depth_size * ((c - 1) + (r - 1) * col_size))
						end
						d := d + 1
					end
					c := c + 1
				end
				r := r + 1
			end
		ensure
			row_set: row_size = big_row
			col_set: col_size = big_col
			depth_set: depth_size = big_depth
		end

	set_spec_zero (br, bc, bd: INTEGER) is
		require
			valid_row: br > 0
			valid_column: bc > 0
			valid_depth: bd > 0
		local
			i: INTEGER
			to_area: TO_SPECIAL [DOUBLE]
		do
			make (br, bc, bd)
			from i := 0 until i >= area.count loop
				area.put (0.0, i)
				i := i + 1
			end
		ensure			
			row_set: row_size = br
			col_set: col_size = bc
			depth_set: depth_size = bd
		end

	make (br, bc, bd: INTEGER) is
		require
			valid_row: br > 0
			valid_column: bc > 0
			valid_depth: bd > 0
		local
			to_area: TO_SPECIAL [DOUBLE]
		do
			row_size := br
			col_size := bc
			depth_size := bd
			!!to_area.make_area (row_size * col_size * depth_size)
			area := to_area.area
			-- !!area.make (row_size * col_size * depth_size)
		end

feature -- size and access

	row_size: INTEGER

	col_size: INTEGER

	depth_size: INTEGER

	item (row, col, dep: INTEGER): DOUBLE is
		require
			valid_row: (row > 0) and (row <= row_size)
			valid_col: (col > 0) and (col <= col_size)
			valid_dep: (dep > 0) and (dep <= depth_size)
		do
			Result := area.item ((dep - 1) + depth_size * ((col - 1) + (row - 1) * col_size))
		end

	put (it: DOUBLE; row, col, dep: INTEGER) is
		require
			valid_row: (row > 0) and (row <= row_size)
			valid_col: (col > 0) and (col <= col_size)
			valid_dep: (dep > 0) and (dep <= depth_size)
		local
			i: INTEGER
		do
			i := (dep - 1) + depth_size * ((col - 1) + (row - 1) * col_size)
			area.put (it, i)
		ensure
			item_there: it = item (row, col, dep)
		end

feature -- conversions

	to_vector (row, col, depth: INTEGER): ARRAY [DOUBLE] is
		require
			row_valid: (row > 0) and (row <= row_size)
			col_valid: (col > 0) and (col <= col_size)
			depth_valid: (depth > 0) and (depth <= depth_size)
		do
			!!Result.make (1, row * col * depth)
			fill_vector (Result, row, col, depth)
		ensure
			result_size: Result.count = row * col * depth
		end

	fill_vector (v: ARRAY [DOUBLE]; row, col, depth: INTEGER) is
		require
			valid_vector: (v /= Void) and then (v.count = row * col * depth)
			row_valid: (row > 0) and (row <= row_size)
			col_valid: (col > 0) and (col <= col_size)
			depth_valid: (depth > 0) and (depth <= depth_size)
		local
			i, r, c, d: INTEGER
			it: DOUBLE
		do
			i := 1
			from r := 1 until r > row loop
				from c := 1 until c > col loop
					from d := 1 until d > depth loop
						it := area.item ((d - 1) + depth_size * ((c - 1) + (r - 1) * col_size))
						-- v.put (item (r, c, d), i)
						v.put (it, i)
						i := i + 1
						d := d + 1
					end
					c := c + 1
				end
				r := r + 1
			end
		end

	fill_and_scale_vector (v: ARRAY [DOUBLE]; factor: DOUBLE; row, col, depth: INTEGER) is
			-- fill vector from data in tensor, and multiply by 
			-- factor when placing data in the vector
		require
			valid_vector: (v /= Void) and then (v.count = row * col * depth)
			row_valid: (row > 0) and (row <= row_size)
			col_valid: (col > 0) and (col <= col_size)
			depth_valid: (depth > 0) and (depth <= depth_size)
		local
			i, r, c, d: INTEGER
			it: DOUBLE
		do
			i := 1
			from r := 1 until r > row loop
				from c := 1 until c > col loop
					from d := 1 until d > depth loop
						it := area.item ((d - 1) + depth_size * ((c - 1) + (r - 1) * col_size))
						-- v.put (item (r, c, d), i)
						v.put (it * factor, i)
						i := i + 1
						d := d + 1
					end
					c := c + 1
				end
				r := r + 1
			end
		end

feature -- operations

	real_element_multiply_by (other: like Current) is
		require
			other.row_size = row_size
			other.col_size = col_size
			other.depth_size = depth_size
		local
			i: INTEGER
			temp: DOUBLE
		do
			from i := 0 until i >= area.count 
			loop
				temp := area.item (i) * other.area.item (i)
				area.put (temp, i)
				i := i + 1
			end
		ensure
			-- result placed in current, other unchanged
		end

	scalar_multiply (x: DOUBLE) is
			-- multiply by a scalar
		local
			i: INTEGER
		do
			from i := 0 until i >= area.count
			loop
				area.put (area.item (i) * x, i)
				i := i + 1
			end
		end

	complex_element_multiply_by (other: like Current) is
		require
			compatible_size: (other.row_size = row_size) and (other.col_size = col_size)
							 and (other.depth_size = depth_size)
			valid_depth: (depth_size \\ 2 = 0)
		local
			r, c, d: INTEGER
			x, y: DOUBLE
			a1, a2, b1, b2: DOUBLE
		do
			from r := 1 until r > row_size loop
				from c := 1 until c > col_size loop
					from d := 1 until d > depth_size loop
						a1 := item (r,c,d)
						a2 := item (r,c,d+1)
						b1 := other.item (r,c,d)
						b2 := other.item(r,c,d+1)
						-- new real part
						-- x := item (r,c,d)*other.item (r,c,d) - item (r,c,d+1)*other.item(r,c,d+1)
						x := a1 * b1 - a2 * b2
						-- new imaginary part
						-- y := item (r,c,d+1)*other.item (r,c,d) + item (r,c,d)*other.item(r,c,d+1)
						y := a2 * b1 + a1 * b2
						-- put (x, r, c, d)
						area.put (x, (d - 1) + depth_size * ((c - 1) + (r - 1) * col_size))
						-- put (y, r, c, d + 1)
						area.put (y, d  + depth_size * ((c - 1) + (r - 1) * col_size))
						d := d + 2
					end
					c := c + 1 
				end
				r := r + 1
			end
		ensure
			-- result placed in current, other unchanged
		end

	conjugate is
			-- replace the tensor by its conjugate (real element 
			-- remains equal, imaginary are * -1)
		require
			depth_even: (depth_size \\ 2) = 0
		local
			r, c, d: INTEGER
			x: DOUBLE
		do
			from r := 1 until r > row_size loop
				from c := 1 until c > col_size loop
					from d := 2 until d > depth_size loop
						x := item (r, c, d)
						put (-x, r, c, d)
						d := d + 2
					end
					c := c + 1
				end
				r := r + 1
			end
		end

	dump is
		local
			i, j , k: INTEGER
			indent: STRING
		do
			print ("FFTW Tensor: ")
			print (row_size)
			print (" x ")
			print  (col_size)
			print (" x ")
			print (depth_size)
			print ("%N")
			indent := " ";
			from i := 1 until (i > depth_size)
			loop
				print ("---Slice: ")
				print (i)
				print ("%N")
				from j := 1 until (j > row_size)
				loop
					print (indent)
					print ("[ ")
					from k := 1 until k > col_size 
					loop
						print (item (j, k, i))
						print (" ")
						k := k + 1
					end
					print (" ]%N")
					j := j + 1
				end
				i := i + 1
			end
		end


feature -- standard features

	copy (other: like Current) is
		do
			row_size := other.row_size
			col_size := other.col_size
			depth_size := other.depth_size
			area := clone (other.area)
		end

	is_equal (other: like Current): BOOLEAN is
		local
			i: INTEGER
		do
			-- check size first
			Result := (row_size = other.row_size)
				and then (col_size = other.col_size)
				and then (depth_size = other.depth_size)
			-- if the same size, compare elements
			if Result then
				from i := 1 until (i >= area.count) or not Result
				loop
					Result := area.item (i) = other.area.item (i)
					i := i + 1
				end
			end
		end

	quasi_is_equal (other: like Current): BOOLEAN is
		require
			other_not_void: other /= Void
		local
			i: INTEGER
			abs: DOUBLE
		do
			Result := (row_size = other.row_size) 
				and then (col_size = other.col_size)
				and then (depth_size = other.depth_size)
			if Result then
				from i := 0 until (i >= area.count) or not Result
				loop
					abs := (area.item (i) - other.area.item (i))
					if abs < 0.0 then 
						abs := - abs
					end
					Result := abs <= 1.0e-15
					i := i + 1
				end
			end
		end


	to_c: ANY is
		do
			Result := area
		end


feature {FFTW_3DTENSOR} -- implementation

	area: SPECIAL [DOUBLE]

end

