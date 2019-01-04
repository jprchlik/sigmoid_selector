make_aia_flare_movies
=====================


Using the same call as get_aia_file_cutout, it creates a movie for every flare with the associated sigmoid at a 45s cadence. In order for this to work, I had to hack aia_mkmovie because aia_prep does not handle cutout files (it blanks the data). I also needed to add a ref_time keyword to aia_mkmovie, which contains anytim format times. The reason for this added keyword is that aia_mkmovie uses files with a specific naming convention. That naming convention is not the same as that delivered by JSOC.