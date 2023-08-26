
This example creates a pts GRIB product, some options are available::

  pts.py --help
  pts.py out1.grib examples/pts/msl_05L_ELSA_2021070300.geo --input points --distance=2.0e5
  pts.py out2.grib 2022062800/pf/*/* --input tc-tracks --filter-wind 62.2
  pts.py out3.grib 2022062800/cf/* --filter-number 1 --filter-basetime "2022-06-28 00" --filter-time 24 72

