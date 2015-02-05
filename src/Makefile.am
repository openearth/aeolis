AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = aeolis
aeolis_SOURCES = constants_module.F90 \
                 sand_mud_bed/modules/message_module.F90 \
                 sand_mud_bed/modules/string_module.F90 \
                 sand_mud_bed/modules/precision_basics.F90 \
                 sand_mud_bed/modules/precision.F90 \
                 sand_mud_bed/modules/bedcomposition_module.F90 \
		 utils_module.F90 \
		 input_module.F90 \
                 output_module.F90 \
		 bed_module.F90 \
                 moist_module.F90 \
		 wind_module.F90 \
		 run_module.F90 \
		 AeoLiS.F90
aeolis_FCFLAGS=-I$(top_srcdir) -I$(srcdir)/sand_mud_bed/include
