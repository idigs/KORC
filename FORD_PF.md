Project: KORC
creation_date:
src_dir: ./src 
output_dir: ./docs
project_github: https://github.com/ORNL-Fusion/KORC  
summary: Kinetic Orbit Runaway electrons Code (KORC)
author: Matt Beidler
author_description: R&amp;D Associate at ORNL
github: https://github.com/mbeidler3
email: beidlermt@ornl.gov
source: true 
graph: true 
search: true 
extensions: f90
print_creation_date: true 
page_dir: User_Documentation
display: public
display: private
display: protected
fixed_length_limit: false

Hi, my name is Matt Beidler.

This is documentation for the Kinetic Orbit Runaway electrons Code
(KORC), developed at Oak Ridge National Laboratory and written
primarily by Leopoldo Carbajal. Mark Cianciosa, Diego
del-Castillo-Negrete and I are currently the primary developers. The
present version of the code follows relativistic electrons in general
electric and magnetic fields under the full Lorentz force, collisions,
and radiation losses.

The present version of the code is compiled and executed on the KNL
nodes of the Cori supercomputer at NERSC.

@Note This documentation is presently under development. Please check
back regularly for updates.

@todo Future plans for development 

1. More physically realistic initial distribution function
2. Implement a guiding-center approximation 
3. Port to multiple computing platforms 
4. Porting to GPU architectures 
5. Coupling with extended-MHD codes to self-consistently evolve 
EM fields and RE distribution functions


@endtodo
