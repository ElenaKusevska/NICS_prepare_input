program NICS_cent
implicit none

!------------------------------------------------------------
! Program to find the centers of rings fo NICS(0) and NICS(1)
! calculation and write some quazi-.gjf file.
!
! Actually, this just reads the input file. The program for
! the NICS calculation is the .cpp file.
!------------------------------------------------------------

integer :: nmols, nrings_total, io, line
character(len=4) :: test_word
integer, allocatable, dimension(:) :: filename_lines, ring_spec_lines
integer, allocatable, dimension(:) :: nrings_per_mol
integer, allocatable, dimension(:,:) :: rings_elements
character(len=20), allocatable, dimension(:) :: filenames
integer :: i, j, k, l, rl


open (unit=3, file='input', status='old', action='read')

! count number of molecules:
nmols = 0
do
   read(3,*,iostat=io) test_word
   if (io/=0) exit
   write(*,*) test_word
   if (test_word == 'mole') then
      nmols = nmols + 1
   end if
end do

allocate (nrings_per_mol(nmols), filename_lines(nmols), filenames(nmols))
do i = 1, nmols
   nrings_per_mol(i) = 0
end do

! count number of rings per molecule and file line where molecule
! names are written:
i = 0 ! molecule number
line = 0 ! file line number
rewind(3)
do
   line = line + 1
   read(3,*,iostat=io) test_word
   if (io/=0) exit
   if (test_word == 'mole') then
      i = i + 1
   end if
   if (test_word == 'file') then
      filename_lines(i) = line
   end if
   if (test_word == 'ring') then
      nrings_per_mol(i) = nrings_per_mol(i) + 1
   end if
end do

nrings_total = sum(nrings_per_mol)

write(*,*) 'filename_lines: ', filename_lines
write(*,*) 'nrings_per_mol: ', nrings_per_mol
write(*,*) 'nrings_total: ', nrings_total 

allocate (ring_spec_lines(nrings_total), rings_elements(nrings_total,6))

! get all filenames:
rewind(3)
i = 0
j = 1
k = 1
rl = 1
do
   i = i + 1
   if (i == filename_lines(j)) then
      read(3,*,iostat=io) test_word, filenames(k)
      k = k + 1
      j = j + 1
   else
      read(3,*,iostat=io) test_word
      if (test_word == 'ring') then
         ring_spec_lines(rl) = i
         rl = rl + 1
      end if
   end if
   if (io/=0) exit
end do

write(*,*) 'ring_spec_lines: ', ring_spec_lines

! get all ring elements:
i = 0
l = 1
rewind(3)
do
   i = i + 1
   if (i == ring_spec_lines(l)) then
      read(3,*,iostat=io) test_word, (rings_elements(l,k), k=1,6)
      l = l + 1
   else
      read(3,*,iostat=io) test_word
      if (io/=0) exit
   end if
   if (l > nrings_total) exit
end do

write(*,*)
write(*,*) filenames
do i = 1, size(rings_elements,1)
   write(*,*) "ring", i, (rings_elements(i,j), j = 1,6)
end do

! sum all elements of nrings_per_mol
! allocate rings_elemens (nrings,6) and initialize to zero
! read input - if ring - read ring elements and increment

close (3)

do i = 1, nmols
   ! open... read filename.xyz
   do j = 1, nrings_per_mol(i)
      ! find center of ring 1
      ! find center of ring 2 ...
      ! write all coordinates to quasi-.gjf file with 0 1 at top
   end do
end do

deallocate(filenames, nrings_per_mol, filename_lines, ring_spec_lines)

end program NICS_cent
