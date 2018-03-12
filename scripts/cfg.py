version_name = "PyFluxPro"
version_number = "V0.1.2"
# V0.1.2 - version produced during resolution of differences
#          between OFQC V2.9.6 and PFP V0.1.1 found by Craig
#          McFarlane at Great Western Woodlands:
#          - auto-complete at L4 was effectively disabled due
#            to a line inserted during debugging being left in
#            place in qcgf.gfalternate_get_correcteddata().
#          - reinstated original line to re-enable auto-complete
# V0.1.1 - cumulative changes up to 13/07/2017
#        - major change to qcutils.CreateSeries()
#          - deprecated "FList=" argument
#          - made all arguments compulsory
#        - QC flag generated for each series immediately
#          prior to qcutils.CreateSeries() based solely
#          on series data mask
# V0.1.0 - copy of OzFluxQC V2.9.6f
