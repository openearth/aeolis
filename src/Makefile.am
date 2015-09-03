AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = aeolis
aeolis_SOURCES = iso_c_utils.F90 \
		 constants_module.F90 \
                 sand_mud_bed/modules/message_module.f90 \
                 sand_mud_bed/modules/string_module.f90 \
                 sand_mud_bed/modules/precision_basics.f90 \
                 sand_mud_bed/modules/precision.f90 \
                 sand_mud_bed/modules/bedcomposition_module.f90 \
		 utils_module.F90 \
		 input_module.F90 \
		 logging.F90 \
                 output_module.F90 \
                 moist_module.F90 \
		 bed_module.F90 \
		 wind_module.F90 \
		 init_module.F90 \
		 run_module.F90 \
		 bmi_module.F90 \
		 aeolis.F90
aeolis_FCFLAGS=-I$(top_srcdir) -I$(srcdir)/sand_mud_bed/include

lib_LTLIBRARIES = libaeolis.la
libaeolis_la_SOURCES = iso_c_utils.F90 \
		 constants_module.F90 \
                 sand_mud_bed/modules/message_module.f90 \
                 sand_mud_bed/modules/string_module.f90 \
                 sand_mud_bed/modules/precision_basics.f90 \
                 sand_mud_bed/modules/precision.f90 \
                 sand_mud_bed/modules/bedcomposition_module.f90 \
		 utils_module.F90 \
		 input_module.F90 \
		 logging.F90 \
                 output_module.F90 \
                 moist_module.F90 \
		 bed_module.F90 \
		 wind_module.F90 \
		 init_module.F90 \
		 run_module.F90 \
		 bmi_module.F90
libaeolis_la_FCFLAGS=-I$(top_srcdir) -I$(srcdir)/sand_mud_bed/include
libaeolis_la_LDFLAGS = -avoid-version -module -shared -dynamiclib -export-dynamic
