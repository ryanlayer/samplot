import unittest
import sys

from samplot import samplot


bam_1 = 'test/data/NA12878_restricted.bam'
bam_2 = 'test/data/NA12889_restricted.bam'
bam_3 = 'test/data/NA12890_restricted.bam'
bams=[bam_1, bam_2, bam_3]

sv_chrm = 'chr4'
sv_start = 115928730
sv_end = 115931875
sv_type = 'DEL'


#{{{ class Test_set_plot_dimensions(unittest.TestCase):
class Test_set_plot_dimensions(unittest.TestCase):
    #{{{ def test_set_plot_dimensions(self):
    def test_set_plot_dimensions(self):

        '''
        def set_plot_dimensions(sv,
                            sv_type,
                            arg_plot_height,
                            arg_plot_width,
                            bams,
                            annotation_files,
                            transcript_file,
                            arg_window,
                            zoom):
        '''
        plot_height = None
        plot_width = None

        annotation_files = None
        transcript_file = None

        zoom = None

        window = None

        sv = [samplot.genome_interval(sv_chrm,sv_start,sv_end)]

        # Test basic function where window is set to be proportional to SV size
        r_plot_height, r_plot_width, r_window, r_ranges = \
            samplot.set_plot_dimensions(sv,
                                        sv_type,
                                        plot_height,
                                        plot_width,
                                        bams,
                                        annotation_files,
                                        transcript_file,
                                        window,
                                        zoom)

        self.assertEqual(r_plot_height, 5)
        self.assertEqual(r_plot_width, 8)
        this_window = int((sv_end - sv_start)/2)
        self.assertEqual( r_window, this_window)
        self.assertEqual( r_ranges[0], 
                          samplot.genome_interval(sv_chrm,
                                                  sv_start - this_window,
                                                  sv_end + this_window))

        # Test to see if zoom is ignored when it is larger than window
        zoom = 10000
        r_plot_height, r_plot_width, r_window, r_ranges = \
            samplot.set_plot_dimensions(sv,
                                        sv_type,
                                        plot_height,
                                        plot_width,
                                        bams,
                                        annotation_files,
                                        transcript_file,
                                        window,
                                        zoom)

        self.assertEqual( r_ranges[0], 
                          samplot.genome_interval(sv_chrm,
                                                  sv_start - this_window,
                                                  sv_end + this_window))


        # Test to see if zoom creates two ranges
        zoom = 100
        r_plot_height, r_plot_width, r_window, r_ranges = \
            samplot.set_plot_dimensions(sv,
                                        sv_type,
                                        plot_height,
                                        plot_width,
                                        bams,
                                        annotation_files,
                                        transcript_file,
                                        window,
                                        zoom)

        self.assertEqual( r_window, zoom)
        self.assertEqual( len(r_ranges), 2)
        self.assertEqual( r_ranges[0], 
                          samplot.genome_interval(sv_chrm,
                                                  sv_start - zoom,
                                                  sv_start + zoom,))
        self.assertEqual( r_ranges[1], 
                          samplot.genome_interval(sv_chrm,
                                                  sv_end - zoom,
                                                  sv_end + zoom) )


        # Test to multiple sv regions
        window = None
        zoom = None
        sv = [samplot.genome_interval(sv_chrm,sv_start,sv_start),
              samplot.genome_interval(sv_chrm,sv_end,sv_end)]
        r_plot_height, r_plot_width, r_window, r_ranges = \
            samplot.set_plot_dimensions(sv,
                                        sv_type,
                                        plot_height,
                                        plot_width,
                                        bams,
                                        annotation_files,
                                        transcript_file,
                                        window,
                                        zoom)

        self.assertEqual( len(r_ranges), 2)
        self.assertEqual( r_ranges[0], 
                          samplot.genome_interval(sv_chrm,
                                                  sv_start-1000,
                                                  sv_start+1000) )
        self.assertEqual( r_ranges[1], 
                          samplot.genome_interval(sv_chrm,
                                                  sv_end-1000,
                                                  sv_end+1000) )
    #}}}

    #{{{def test_get_read_data(self):
    def test_get_read_data(self):
        '''
        read_data,max_coverage = get_read_data(ranges,
                                               options.bams,
                                               options.reference,
                                               options.min_mqual,
                                               options.coverage_only,
                                               options.long_read,
                                               options.same_yaxis_scales,
                                               options.max_depth,
                                               options.z)
        '''

        plot_height = None
        plot_width = None

        annotation_files = None
        transcript_file = None

        zoom = None

        window = None

        sv = [samplot.genome_interval(sv_chrm,sv_start,sv_end)]

        # Test basic function where window is set to be proportional to SV size
        r_plot_height, r_plot_width, r_window, r_ranges = \
            samplot.set_plot_dimensions(sv,
                                        sv_type,
                                        plot_height,
                                        plot_width,
                                        bams,
                                        annotation_files,
                                        transcript_file,
                                        window,
                                        zoom)

        reference = None
        min_mqual = None
        coverage_only = None
        long_read = 1000
        long_even_size = 100
        same_yaxis_scales = None
        max_depth = 100
        z = 4

        read_data,max_coverage = samplot.get_read_data(r_ranges,
                                                       bams,
                                                       reference,
                                                       min_mqual,
                                                       coverage_only,
                                                       long_read,
                                                       long_even_size,
                                                       same_yaxis_scales,
                                                       max_depth,
                                                       z)
    #}}}
#}}}

#{{{ class Test_genome_interval(unittest.TestCase):
class Test_genome_interval(unittest.TestCase):
    #{{{ def test_init(self):
    def test_init(self):
        gi = samplot.genome_interval('chr1', 1, 1000)
        self.assertEqual(gi.chrm, 'chr1')
        self.assertEqual(gi.start, 1)
        self.assertEqual(gi.end, 1000)
    #}}}

    #{{{ def test_init(self):
    def test_intersect(self):
        gi = samplot.genome_interval('chr8', 500, 1000)

        self.assertEqual(-1, 
                         gi.intersect(samplot.genome_interval('chr7',
                                                               500,
                                                               1000)))

        self.assertEqual(1,
                         gi.intersect(samplot.genome_interval('chr9',
                                                              500,
                                                              1000)))

        self.assertEqual(-1,
                         gi.intersect(samplot.genome_interval('chr8',
                                                              100,
                                                              499)))

        self.assertEqual(1,
                         gi.intersect(samplot.genome_interval('chr8',
                                                              1001,
                                                              2000)))

        self.assertEqual(0,
                         gi.intersect(samplot.genome_interval('chr8',
                                                              1,
                                                              500)))
        self.assertEqual(0,
                         gi.intersect(samplot.genome_interval('chr8',
                                                              500,
                                                              501)))
        self.assertEqual(0,
                         gi.intersect(samplot.genome_interval('chr8',
                                                              1000,
                                                              2000)))

    #}}}

    #{{{ def test_get_range_hit(self):
    def test_get_range_hit(self):
        gi_0 = samplot.genome_interval('chr8', 500, 1000)
        ranges = [gi_0]

        self.assertEqual(0, samplot.get_range_hit(ranges, 'chr8', 500))


        gi_1 = samplot.genome_interval('chr8', 2000, 3000)
        ranges = [gi_0, gi_1]
        self.assertEqual(0, samplot.get_range_hit(ranges, 'chr8', 500))
        self.assertEqual(1, samplot.get_range_hit(ranges, 'chr8', 2500))

        self.assertEqual(None, samplot.get_range_hit(ranges, 'chr7', 2500))
        self.assertEqual(None, samplot.get_range_hit(ranges, 'chr8', 100))
        self.assertEqual(None, samplot.get_range_hit(ranges, 'chr8', 10000))
    #}}}

    #{{{ def test_map_genome_point_to_range_points(self):
    def test_map_genome_point_to_range_points(self):
        gi_0 = samplot.genome_interval('chr8', 100, 200)
        ranges = [gi_0]

        self.assertEqual(None,
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  10))
        self.assertEqual(0.0, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  100))
        self.assertEqual(0.25, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  125))
        self.assertEqual(0.5, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  150))
        self.assertEqual(0.75, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  175))
        self.assertEqual(1.0, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  200))
        self.assertEqual(None,
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  201))

        gi_1 = samplot.genome_interval('chr8', 300, 400)
        ranges = [gi_0, gi_1]

        self.assertEqual(None,
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  10))
        self.assertEqual(0.0, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  100))
        self.assertEqual(0.25/2, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  125))
        self.assertEqual(0.5/2, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  150))
        self.assertEqual(0.75/2, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  175))
        self.assertEqual(1.0/2, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  200))
        self.assertEqual(None,
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  201))
        self.assertEqual(0.5, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  300))
        self.assertEqual(0.5+0.25/2, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  325))
        self.assertEqual(0.5+0.5/2, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  350))
        self.assertEqual(0.5+0.75/2, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  375))
        self.assertEqual(1.0, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  400))
        gi_0 = samplot.genome_interval('chr8', 100, 200)
        gi_1 = samplot.genome_interval('chr9', 300, 400)
        ranges = [gi_0, gi_1]

        self.assertEqual(None,
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  10))
        self.assertEqual(0.0, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  100))
        self.assertEqual(0.25/2, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  125))
        self.assertEqual(0.5/2, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  150))
        self.assertEqual(0.75/2, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  175))
        self.assertEqual(1.0/2, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  200))
        self.assertEqual(None,
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr8',
                                                                  201))
        self.assertEqual(0.5, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr9',
                                                                  300))
        self.assertEqual(0.5+0.25/2, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr9',
                                                                  325))
        self.assertEqual(0.5+0.5/2, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr9',
                                                                  350))
        self.assertEqual(0.5+0.75/2, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr9',
                                                                  375))
        self.assertEqual(1.0, 
                         samplot.map_genome_point_to_range_points(ranges,
                                                                  'chr9',
                                                                  400))
    #}}}

#}}}

#{{{ class Test_long_read_plan(unittest.TestCase):
class Test_long_read_plan(unittest.TestCase):
    #{{{ def test_init(self):
    def test_add_align_step(self):
        alignment = samplot.Alignment('chr8', 100, 500, True, 0)

        # both are in the same range
        gi_0 = samplot.genome_interval('chr8', 100, 1000)
        ranges = [gi_0]
        steps = []

        samplot.add_align_step(alignment, steps, ranges)

        self.assertEqual(1, len(steps))
        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(100, steps[0].start_pos.start)
        self.assertEqual(100, steps[0].start_pos.end)
        self.assertEqual('chr8', steps[0].end_pos.chrm)
        self.assertEqual(500, steps[0].end_pos.start)
        self.assertEqual(500, steps[0].end_pos.end)
        self.assertEqual('Align', steps[0].info['TYPE'])

        # in different ranges
        gi_0 = samplot.genome_interval('chr8', 100, 200)
        gi_1 = samplot.genome_interval('chr8', 300, 1000)
        ranges = [gi_0, gi_1]
        steps = []

        samplot.add_align_step(alignment, steps, ranges)

        self.assertEqual(2, len(steps))
        #start
        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(100, steps[0].start_pos.start)
        self.assertEqual(100, steps[0].start_pos.end)
        #end
        self.assertEqual('chr8', steps[0].end_pos.chrm)
        self.assertEqual(200, steps[0].end_pos.start)
        self.assertEqual(200, steps[0].end_pos.end)
        #event
        self.assertEqual('Align', steps[0].info['TYPE'])

        #start
        self.assertEqual('chr8', steps[1].start_pos.chrm)
        self.assertEqual(300, steps[1].start_pos.start)
        self.assertEqual(300, steps[1].start_pos.end)
        #end
        self.assertEqual('chr8', steps[1].end_pos.chrm)
        self.assertEqual(500, steps[1].end_pos.start)
        self.assertEqual(500, steps[1].end_pos.end)
        #event
        self.assertEqual('Align', steps[1].info['TYPE'])

        # start is not in range, use end hit
        gi_0 = samplot.genome_interval('chr8', 10, 20)
        gi_1 = samplot.genome_interval('chr8', 300, 1000)
        ranges = [gi_0, gi_1]
        steps = []

        samplot.add_align_step(alignment, steps, ranges)

        self.assertEqual(1, len(steps))
        #start
        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(300, steps[0].start_pos.start)
        self.assertEqual(300, steps[0].start_pos.end)
        #end
        self.assertEqual('chr8', steps[0].end_pos.chrm)
        self.assertEqual(500, steps[0].end_pos.start)
        self.assertEqual(500, steps[0].end_pos.end)
        #event
        self.assertEqual('Align', steps[0].info['TYPE'])

        # end is not in range, use start hit
        gi_0 = samplot.genome_interval('chr8', 100, 200)
        gi_1 = samplot.genome_interval('chr8', 3000, 4000)
        ranges = [gi_0, gi_1]
        steps = []

        samplot.add_align_step(alignment, steps, ranges)

        #start
        self.assertEqual(1, len(steps))
        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(100, steps[0].start_pos.start)
        self.assertEqual(100, steps[0].start_pos.end)
        #end
        self.assertEqual('chr8', steps[0].end_pos.chrm)
        self.assertEqual(200, steps[0].end_pos.end)
        self.assertEqual(200, steps[0].end_pos.start)
        #event
        self.assertEqual('Align', steps[0].info['TYPE'])

        # neither end is in range, add nothing
        gi_0 = samplot.genome_interval('chr8', 10, 20)
        gi_1 = samplot.genome_interval('chr8', 3000, 4000)
        ranges = [gi_0, gi_1]
        steps = []

        samplot.add_align_step(alignment, steps, ranges)

        self.assertEqual(0, len(steps))
    #}}}
     
    #{{{def test_get_alignments_from_cigar(self):
    def test_get_alignments_from_cigar(self):
        '''
        alignments = get_alignments_from_cigar(
                bam_file.get_reference_name(read.reference_id),
                read.pos,
                not read.is_reverse,
                read.cigartuples)
        '''
        CIGAR_MAP = { 'M' : 0,
                      'I' : 1,
                      'D' : 2,
                      'N' : 3,
                      'S' : 4,
                      'H' : 5,
                      'P' : 6,
                      '=' : 7,
                      'X' : 8,
                      'B' : 9 }

        cigar = [(CIGAR_MAP['M'], 100),
                 (CIGAR_MAP['D'], 100),
                 (CIGAR_MAP['M'], 100)]
        alignments = samplot.get_alignments_from_cigar('chr8',
                                                       100,
                                                       True,
                                                       cigar)
        self.assertEqual(2,len(alignments))

        self.assertEqual('chr8', alignments[0].pos.chrm)
        self.assertEqual(100, alignments[0].pos.start)
        self.assertEqual(200, alignments[0].pos.end)
        self.assertEqual(True, alignments[0].strand)
        self.assertEqual(0, alignments[0].query_position)

        self.assertEqual('chr8', alignments[1].pos.chrm)
        self.assertEqual(300, alignments[1].pos.start)
        self.assertEqual(400, alignments[1].pos.end)
        self.assertEqual(True, alignments[1].strand)
        self.assertEqual(100, alignments[1].query_position)
    #}}}

    #{{{def test_get_long_read_plan(self):
    def test_get_long_read_plan(self):
        gi_0 = samplot.genome_interval('chr8', 100, 250)
        gi_1 = samplot.genome_interval('chr8', 300, 400)
        ranges = [gi_0, gi_1]
        long_reads = {}
        read_name = 'Test'
        alignments = [samplot.Alignment('chr8', 100, 200, True, 0)]
        long_reads[read_name] = [ samplot.LongRead(alignments) ]


        max_gap, steps = samplot.get_long_read_plan(read_name,
                                                    long_reads,
                                                    ranges)

        self.assertEqual(0, max_gap)
        self.assertEqual(1, len(steps))
        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(100, steps[0].start_pos.start)
        self.assertEqual(100, steps[0].start_pos.end)
        self.assertEqual('chr8', steps[0].end_pos.chrm)
        self.assertEqual(200, steps[0].end_pos.start)
        self.assertEqual(200, steps[0].end_pos.end)
        self.assertEqual('LONGREAD', steps[0].event)
        self.assertEqual('Align', steps[0].info['TYPE'])

        alignments = [samplot.Alignment('chr8', 100, 299, True, 0)]
        long_reads[read_name] = [ samplot.LongRead(alignments) ]
        max_gap, steps = samplot.get_long_read_plan(read_name,
                                                    long_reads,
                                                    ranges)


        self.assertEqual(0, max_gap)
        self.assertEqual(1, len(steps))
        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(100, steps[0].start_pos.start)
        self.assertEqual(100, steps[0].start_pos.end)
        self.assertEqual('chr8', steps[0].end_pos.chrm)
        self.assertEqual(250, steps[0].end_pos.start)
        self.assertEqual(250, steps[0].end_pos.end)
        self.assertEqual('Align', steps[0].info['TYPE'])


        alignments = [samplot.Alignment('chr8', 100, 350, True, 0)]
        long_reads[read_name] = [ samplot.LongRead(alignments) ]
        max_gap, steps = samplot.get_long_read_plan(read_name,
                                                    long_reads,
                                                    ranges)

        self.assertEqual(0, max_gap)
        self.assertEqual(2, len(steps))
        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(100, steps[0].start_pos.start)
        self.assertEqual(100, steps[0].start_pos.end)
        self.assertEqual('chr8', steps[0].end_pos.chrm)
        self.assertEqual(250, steps[0].end_pos.start)
        self.assertEqual(250, steps[0].end_pos.end)
        self.assertEqual('Align', steps[0].info['TYPE'])

        self.assertEqual('chr8', steps[1].start_pos.chrm)
        self.assertEqual(300, steps[1].start_pos.start)
        self.assertEqual(300, steps[1].start_pos.end)
        self.assertEqual('chr8', steps[1].end_pos.chrm)
        self.assertEqual(350, steps[1].end_pos.start)
        self.assertEqual(350, steps[1].end_pos.end)
        self.assertEqual('Align', steps[1].info['TYPE'])

        alignments = [samplot.Alignment('chr8', 100, 250, True, 0),
                      samplot.Alignment('chr8', 300, 350, True, 150)]
        long_reads[read_name] = [ samplot.LongRead(alignments) ]
        max_gap, steps = samplot.get_long_read_plan(read_name,
                                                    long_reads,
                                                    ranges)

        self.assertEqual(50, max_gap)
        self.assertEqual(3, len(steps))

        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(100, steps[0].start_pos.start)
        self.assertEqual(100, steps[0].start_pos.end)
        self.assertEqual('chr8', steps[0].end_pos.chrm)
        self.assertEqual(250, steps[0].end_pos.start)
        self.assertEqual(250, steps[0].end_pos.end)
        self.assertEqual('Align', steps[0].info['TYPE'])

        self.assertEqual('chr8', steps[1].start_pos.chrm)
        self.assertEqual(250, steps[1].start_pos.start)
        self.assertEqual(250, steps[1].start_pos.end)
        self.assertEqual('chr8', steps[1].end_pos.chrm)
        self.assertEqual(300, steps[1].end_pos.start)
        self.assertEqual(300, steps[1].end_pos.end)
        self.assertEqual('Deletion', steps[1].info['TYPE'])

        self.assertEqual('chr8', steps[2].start_pos.chrm)
        self.assertEqual(300, steps[2].start_pos.start)
        self.assertEqual(300, steps[2].start_pos.end)
        self.assertEqual('chr8', steps[2].end_pos.chrm)
        self.assertEqual(350, steps[2].end_pos.start)
        self.assertEqual(350, steps[2].end_pos.end)
        self.assertEqual('Align', steps[2].info['TYPE'])

        gi_0 = samplot.genome_interval('chr8', 100, 250)
        gi_1 = samplot.genome_interval('chr9', 300, 400)
        ranges = [gi_0, gi_1]
 
        alignments = [samplot.Alignment('chr8', 100, 250, True, 0),
                      samplot.Alignment('chr9', 300, 350, True, 150)]
        long_reads[read_name] = [ samplot.LongRead(alignments) ]
        max_gap, steps = samplot.get_long_read_plan(read_name,
                                                    long_reads,
                                                    ranges)

        self.assertEqual(5000, max_gap)
        self.assertEqual(3, len(steps))

        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(100, steps[0].start_pos.start)
        self.assertEqual(100, steps[0].start_pos.end)
        self.assertEqual('chr8', steps[0].end_pos.chrm)
        self.assertEqual(250, steps[0].end_pos.start)
        self.assertEqual(250, steps[0].end_pos.end)
        self.assertEqual('Align', steps[0].info['TYPE'])

        self.assertEqual('chr8', steps[1].start_pos.chrm)
        self.assertEqual(250, steps[1].start_pos.start)
        self.assertEqual(250, steps[1].start_pos.end)
        self.assertEqual('chr9', steps[1].end_pos.chrm)
        self.assertEqual(300, steps[1].end_pos.start)
        self.assertEqual(300, steps[1].end_pos.end)
        self.assertEqual('InterChrm', steps[1].info['TYPE'])

        self.assertEqual('chr9', steps[2].start_pos.chrm)
        self.assertEqual(300, steps[2].start_pos.start)
        self.assertEqual(300, steps[2].start_pos.end)
        self.assertEqual('chr9', steps[2].end_pos.chrm)
        self.assertEqual(350, steps[2].end_pos.start)
        self.assertEqual(350, steps[2].end_pos.end)
        self.assertEqual('Align', steps[2].info['TYPE'])
    #}}}
#}}}

#{{{class Test_annotation_plan(unittest.TestCase):
class Test_annotation_plan(unittest.TestCase):
    #{{{def test_get_alignments_from_cigar(self):
    def test_get_alignments_from_cigar(self):

        gi_1 = samplot.genome_interval('chr8', 100, 200)
        gi_2 = samplot.genome_interval('chr8', 300, 400)
        ranges = [gi_1, gi_2]

        i = samplot.genome_interval('chr8', 110, 120)
        s, e = samplot.get_interval_range_plan_start_end(ranges, i)

        self.assertEqual('chr8',s.chrm)
        self.assertEqual(110,s.start)
        self.assertEqual(110,s.end)
        self.assertEqual('chr8',e.chrm)
        self.assertEqual(120,e.start)
        self.assertEqual(120,e.end)


        i = samplot.genome_interval('chr8', 110, 220)
        s, e = samplot.get_interval_range_plan_start_end(ranges, i)

        self.assertEqual('chr8',s.chrm)
        self.assertEqual(110,s.start)
        self.assertEqual(110,s.end)
        self.assertEqual('chr8',e.chrm)
        self.assertEqual(200,e.start)
        self.assertEqual(200,e.end)


        i = samplot.genome_interval('chr8', 220, 320)
        s, e = samplot.get_interval_range_plan_start_end(ranges, i)

        self.assertEqual('chr8',s.chrm)
        self.assertEqual(300,s.start)
        self.assertEqual(300,s.end)
        self.assertEqual('chr8',e.chrm)
        self.assertEqual(320,e.start)
        self.assertEqual(320,e.end)


        i = samplot.genome_interval('chr8', 120, 320)
        s, e = samplot.get_interval_range_plan_start_end(ranges, i)

        self.assertEqual('chr8',s.chrm)
        self.assertEqual(120,s.start)
        self.assertEqual(120,s.end)
        self.assertEqual('chr8',e.chrm)
        self.assertEqual(320,e.start)
        self.assertEqual(320,e.end)

        i = samplot.genome_interval('chr8', 320, 520)
        s, e = samplot.get_interval_range_plan_start_end(ranges, i)

        self.assertEqual('chr8',s.chrm)
        self.assertEqual(320,s.start)
        self.assertEqual(320,s.end)
        self.assertEqual('chr8',e.chrm)
        self.assertEqual(400,e.start)
        self.assertEqual(400,e.end)


        i = samplot.genome_interval('chr8', 30, 50)
        s, e = samplot.get_interval_range_plan_start_end(ranges, i)

        self.assertEqual(None, s)
        self.assertEqual(None, e)

        i = samplot.genome_interval('chr8', 3000, 5000)
        s, e = samplot.get_interval_range_plan_start_end(ranges, i)

        self.assertEqual(None, s)
        self.assertEqual(None, e)
    #}}}
#}}}

#{{{class Test_splits(unittest.TestCase):
class Test_splits(unittest.TestCase):
    #{{{def test_get_split_plan(self):
    def test_get_split_plan(self):

        splits = {}
        hp = 0
        splits[hp] = {}
        read_name_1 = 'Test1'

        ranges = [samplot.genome_interval('chr8', 100, 200),
                  samplot.genome_interval('chr8', 600, 800) ]

        #both in same ragne
        #Deletion
        splits[hp][read_name_1] = [\
                samplot.SplitRead('chr8', 100, 150, True, 0, False, False),
                samplot.SplitRead('chr8', 170, 180, True, 50, False, False)]

        plan = samplot.get_split_plan(ranges, splits[hp][read_name_1])

        max_gap, steps = plan

        self.assertEqual(20, max_gap)
        self.assertEqual(1, len(steps))
        self.assertEqual('SPLITREAD', steps[0].event)
        self.assertEqual('Deletion', steps[0].info['TYPE'])
        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(150, steps[0].start_pos.start)
        self.assertEqual(150, steps[0].start_pos.end)
        self.assertEqual('chr8', steps[0].end_pos.chrm)
        self.assertEqual(170, steps[0].end_pos.start)
        self.assertEqual(170, steps[0].end_pos.end)

        #Duplication
        splits[hp][read_name_1] = [\
                samplot.SplitRead('chr8', 100, 150, True, 0, False, False),
                samplot.SplitRead('chr8', 130, 180, True, 50, False, False)]

        plan = samplot.get_split_plan(ranges, splits[hp][read_name_1])

        max_gap, steps = plan

        self.assertEqual(20, max_gap)
        self.assertEqual(1, len(steps))
        self.assertEqual('SPLITREAD', steps[0].event)
        self.assertEqual('Duplication', steps[0].info['TYPE'])
        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(150, steps[0].start_pos.start)
        self.assertEqual(150, steps[0].start_pos.end)
        self.assertEqual('chr8', steps[0].end_pos.chrm)
        self.assertEqual(130, steps[0].end_pos.start)
        self.assertEqual(130, steps[0].end_pos.end)

        #Inversion
        splits[hp][read_name_1] = [\
                samplot.SplitRead('chr8', 100, 150, True, 0, False, False),
                samplot.SplitRead('chr8', 151, 180, False, 50, False, False)]

        plan = samplot.get_split_plan(ranges, splits[hp][read_name_1])

        max_gap, steps = plan

        self.assertEqual(30, max_gap)
        self.assertEqual(1, len(steps))
        self.assertEqual('SPLITREAD', steps[0].event)
        self.assertEqual('Inversion', steps[0].info['TYPE'])
        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(150, steps[0].start_pos.start)
        self.assertEqual(150, steps[0].start_pos.end)
        self.assertEqual('chr8', steps[0].end_pos.chrm)
        self.assertEqual(151, steps[0].end_pos.start)
        self.assertEqual(151, steps[0].end_pos.end)


        #both in same ragne
        splits[hp][read_name_1] = [\
                samplot.SplitRead('chr8', 100, 150, True, 0, False, False)]

        plan = samplot.get_split_plan(ranges, splits[hp][read_name_1])
        self.assertEqual(None, plan)

        #both in same ragne
        splits[hp][read_name_1] = [\
                samplot.SplitRead('chr8', 550, 650, True, 0, False, False),
                samplot.SplitRead('chr8', 700, 750, True, 50, False, False)]

        plan = samplot.get_split_plan(ranges, splits[hp][read_name_1])

        max_gap, steps = plan

        self.assertEqual(50, max_gap)
        self.assertEqual(1, len(steps))
        self.assertEqual('SPLITREAD', steps[0].event)
        self.assertEqual('Deletion', steps[0].info['TYPE'])
        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(650, steps[0].start_pos.start)
        self.assertEqual(650, steps[0].start_pos.end)
        self.assertEqual('chr8', steps[0].end_pos.chrm)
        self.assertEqual(700, steps[0].end_pos.start)
        self.assertEqual(700, steps[0].end_pos.end)



        #both in same ragne
        splits[hp][read_name_1] = [\
                samplot.SplitRead('chr8', 150, 175, True, 0, False, False),
                samplot.SplitRead('chr8', 650, 675, True, 50, False, False)]

        plan = samplot.get_split_plan(ranges, splits[hp][read_name_1])

        max_gap, steps = plan

        self.assertEqual(475, max_gap)
        self.assertEqual(1, len(steps))
        self.assertEqual('SPLITREAD', steps[0].event)
        self.assertEqual('Deletion', steps[0].info['TYPE'])
        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(175, steps[0].start_pos.start)
        self.assertEqual(175, steps[0].start_pos.end)
        self.assertEqual('chr8', steps[0].end_pos.chrm)
        self.assertEqual(650, steps[0].end_pos.start)
        self.assertEqual(650, steps[0].end_pos.end)


        #inter chrom
        ranges = [samplot.genome_interval('chr8', 100, 200),
                  samplot.genome_interval('chr9', 600, 800) ]

        splits[hp][read_name_1] = [\
                samplot.SplitRead('chr8', 150, 175, True, 0, False, False),
                samplot.SplitRead('chr9', 650, 675, True, 50, False, False)]

        plan = samplot.get_split_plan(ranges, splits[hp][read_name_1])

        max_gap, steps = plan

        self.assertEqual(samplot.INTERCHROM_YAXIS, max_gap)
        self.assertEqual(1, len(steps))
        self.assertEqual('SPLITREAD', steps[0].event)
        self.assertEqual('InterChrm', steps[0].info['TYPE'])
        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(175, steps[0].start_pos.start)
        self.assertEqual(175, steps[0].start_pos.end)
        self.assertEqual('chr9', steps[0].end_pos.chrm)
        self.assertEqual(650, steps[0].end_pos.start)
        self.assertEqual(650, steps[0].end_pos.end)


        splits[hp][read_name_1] = [\
                samplot.SplitRead('chr8', 150, 175, True, 0, False, False),
                samplot.SplitRead('chr9', 650, 675, False, 50, False, False)]

        plan = samplot.get_split_plan(ranges, splits[hp][read_name_1])

        max_gap, steps = plan

        self.assertEqual(samplot.INTERCHROM_YAXIS, max_gap)
        self.assertEqual(1, len(steps))
        self.assertEqual('SPLITREAD', steps[0].event)
        self.assertEqual('InterChrmInversion', steps[0].info['TYPE'])
        self.assertEqual('chr8', steps[0].start_pos.chrm)
        self.assertEqual(175, steps[0].start_pos.start)
        self.assertEqual(175, steps[0].start_pos.end)
        self.assertEqual('chr9', steps[0].end_pos.chrm)
        self.assertEqual(650, steps[0].end_pos.start)
        self.assertEqual(650, steps[0].end_pos.end)
    #}}}

    #{{{def test_get_splits_plan(self):
    def test_get_splits_plan(self):

        splits = {}
        hp = 0
        splits[hp] = {}

        ranges = [samplot.genome_interval('chr8', 100, 200),
                  samplot.genome_interval('chr9', 600, 800) ]

        #Deletion
        splits[hp]['del'] = [\
                samplot.SplitRead('chr8', 100, 150, True, 0, False, False),
                samplot.SplitRead('chr8', 170, 180, True, 50, False, False)]

        #Duplication
        splits[hp]['dup'] = [\
                samplot.SplitRead('chr8', 100, 150, True, 0, False, False),
                samplot.SplitRead('chr8', 130, 180, True, 50, False, False)]

        #Inversion
        splits[hp]['inv'] = [\
                samplot.SplitRead('chr8', 100, 150, True, 0, False, False),
                samplot.SplitRead('chr8', 151, 180, False, 50, False, False)]

        #Bad split
        splits[hp]['bad'] = [\
                samplot.SplitRead('chr8', 100, 150, True, 0, False, False)]


        #Interchm
        splits[hp]['interchm'] = [\
                samplot.SplitRead('chr8', 150, 175, True, 0, False, False),
                samplot.SplitRead('chr9', 650, 675, True, 50, False, False)]

        #InterchmInv
        splits[hp]['interchminv'] = [\
                samplot.SplitRead('chr8', 150, 175, True, 0, False, False),
                samplot.SplitRead('chr9', 650, 675, False, 50, False, False)]

        plan = samplot.get_splits_plan(ranges, splits[hp])

        max_gap, steps = plan

        self.assertEqual(samplot.INTERCHROM_YAXIS, max_gap)
        self.assertEqual(5, len(steps))
    #}}}


   #{{{def test_get_split_insert_size(self):
    def test_get_split_insert_size(self):

        splits = {}
        hp = 0
        splits[hp] = {}
        read_name_1 = 'Test1'

        #both in same ragne
        splits[hp][read_name_1] = [\
                samplot.SplitRead('chr8', 100, 150, True, 0, False, False),
                samplot.SplitRead('chr8', 170, 180, True, 50, False, False)]

        read_name_2 = 'Test2'

        splits[hp][read_name_2] = [\
                samplot.SplitRead('chr8', 100, 150, True, 0, False, False),
                samplot.SplitRead('chr8', 170, 180, True, 50, False, False)]

        ranges = [samplot.genome_interval('chr8', 100, 200),
                  samplot.genome_interval('chr8', 600, 800) ]

        split_insert_sizes = samplot.get_splits_insert_sizes(ranges, splits)

        self.assertEqual(2, len(split_insert_sizes))
        self.assertEqual(20, split_insert_sizes[0])
        self.assertEqual(20, split_insert_sizes[1])

        #one starting in range ends out of range
        splits[hp][read_name_1] = [\
                samplot.SplitRead('chr8', 100, 350, True, 0, False, False),
                samplot.SplitRead('chr8', 650, 700, True, 250, False, False)]

        split_insert_sizes = samplot.get_splits_insert_sizes(ranges, splits)

        self.assertEqual(2, len(split_insert_sizes))
        self.assertEqual(300, split_insert_sizes[0])
        self.assertEqual(20, split_insert_sizes[1])

        #one out of range
        splits[hp][read_name_1] = [\
                samplot.SplitRead('chr8', 10, 35, True, 0, False, False),
                samplot.SplitRead('chr8', 650, 700, True, 25, False, False)]

        split_insert_sizes = samplot.get_splits_insert_sizes(ranges, splits)

        self.assertEqual(1, len(split_insert_sizes))
        self.assertEqual(20, split_insert_sizes[0])

        #DUP
        splits[hp][read_name_1] = [\
                samplot.SplitRead('chr8', 125, 150, True, 0, False, False),
                samplot.SplitRead('chr8', 130, 155, True, 25, False, False)]
        split_insert_sizes = samplot.get_splits_insert_sizes(ranges, splits)
        self.assertEqual(2, len(split_insert_sizes))
        self.assertEqual(20, split_insert_sizes[0])
        self.assertEqual(20, split_insert_sizes[1])

        #INV
        splits[hp][read_name_1] = [\
                samplot.SplitRead('chr8', 125, 150, True, 0, False, False),
                samplot.SplitRead('chr8', 151, 175, False, 25, False, False)]
        split_insert_sizes = samplot.get_splits_insert_sizes(ranges, splits)
        self.assertEqual(2, len(split_insert_sizes))
        self.assertEqual(25, split_insert_sizes[0])
        self.assertEqual(20, split_insert_sizes[1])

        #interchrm
        ranges = [samplot.genome_interval('chr8', 100, 200),
                  samplot.genome_interval('chr9', 600, 800) ]
        splits[hp][read_name_1] = [\
                samplot.SplitRead('chr8', 125, 150, True, 0, False, False),
                samplot.SplitRead('chr9', 650, 675, True, 25, False, False)]
        split_insert_sizes = samplot.get_splits_insert_sizes(ranges, splits)
        self.assertEqual(2, len(split_insert_sizes))
        self.assertEqual(samplot.INTERCHROM_YAXIS, split_insert_sizes[0])
        self.assertEqual(20, split_insert_sizes[1])


    #}}}
#}}}

#{{{ class Test_pairs(unittest.TestCase):
class Test_pairs(unittest.TestCase):
    #{{{ def test_get_pair_insert_size(self):
    def test_get_pair_insert_size(self):

        ranges = [samplot.genome_interval('chr8', 100, 200),
                  samplot.genome_interval('chr8', 600, 800) ]

        pairs = {}
        hp = 0
        pairs[hp] = {}
        read_name_1 = 'Test1'

        #both in same ragne
        pairs[hp][read_name_1] = [\
                samplot.PairedEnd('chr8', 100, 150, True, False, False),
                samplot.PairedEnd('chr8', 170, 180, False, False, False)]

        read_name_2 = 'Test2'

        pairs[hp][read_name_2] = [\
                samplot.PairedEnd('chr8', 100, 150, True, False, False),
                samplot.PairedEnd('chr8', 170, 180, False, False, False)]
        pair_insert_sizes = samplot.get_pairs_insert_sizes(ranges, pairs)

        self.assertEqual(2, len(pair_insert_sizes))
        self.assertEqual(80, pair_insert_sizes[0])
        self.assertEqual(80, pair_insert_sizes[1])

        #one starting in range ends out of range
        pairs[hp][read_name_1] = [\
                samplot.PairedEnd('chr8', 100, 150, True, False, False),
                samplot.PairedEnd('chr8', 190, 240, False, False, False)]

        pair_insert_sizes = samplot.get_pairs_insert_sizes(ranges, pairs)

        self.assertEqual(2, len(pair_insert_sizes))
        self.assertEqual(140, pair_insert_sizes[0])
        self.assertEqual(80, pair_insert_sizes[1])


        #one out of range
        pairs[hp][read_name_1] = [\
                samplot.PairedEnd('chr9', 100, 150, True, False, False),
                samplot.PairedEnd('chr8', 190, 240, False, False, False)]

        pair_insert_sizes = samplot.get_pairs_insert_sizes(ranges, pairs)

        self.assertEqual(1, len(pair_insert_sizes))
        self.assertEqual(80, pair_insert_sizes[0])


        #DUP
        pairs[hp][read_name_1] = [\
                samplot.PairedEnd('chr8', 125, 150, True, False, False),
                samplot.PairedEnd('chr8', 175, 200, False, False, False)]

        pair_insert_sizes = samplot.get_pairs_insert_sizes(ranges, pairs)

        self.assertEqual(2, len(pair_insert_sizes))
        self.assertEqual(75, pair_insert_sizes[0])
        self.assertEqual(80, pair_insert_sizes[1])


        #INV
        pairs[hp][read_name_1] = [\
                samplot.PairedEnd('chr8', 125, 150, True, False, False),
                samplot.PairedEnd('chr8', 175, 200, True, False, False)]

        pair_insert_sizes = samplot.get_pairs_insert_sizes(ranges, pairs)

        self.assertEqual(2, len(pair_insert_sizes))
        self.assertEqual(75, pair_insert_sizes[0])
        self.assertEqual(80, pair_insert_sizes[1])

        #interchrm
        ranges = [samplot.genome_interval('chr8', 100, 200),
                  samplot.genome_interval('chr9', 600, 800) ]

        pairs[hp][read_name_1] = [\
                samplot.PairedEnd('chr8', 125, 150, True, False, False),
                samplot.PairedEnd('chr9', 675, 700, True, False, False)]

        pair_insert_sizes = samplot.get_pairs_insert_sizes(ranges, pairs)

        self.assertEqual(2, len(pair_insert_sizes))
        self.assertEqual(samplot.INTERCHROM_YAXIS, pair_insert_sizes[0])
        self.assertEqual(80, pair_insert_sizes[1])

    #}}}

    #{{{ def test_get_pair_plan(self):
    def test_get_pair_plan(self):

        ranges = [samplot.genome_interval('chr8', 100, 200),
                  samplot.genome_interval('chr8', 600, 800) ]

        pairs = {}
        hp = 0
        pairs[hp] = {}
        read_name_1 = 'Test1'

        #both in same ragne
        pairs[hp][read_name_1] = [\
                samplot.PairedEnd('chr8', 100, 150, False, False, False),
                samplot.PairedEnd('chr8', 170, 180, True, False, False)]

        read_name_2 = 'Test2'

        pairs[hp][read_name_2] = [\
                samplot.PairedEnd('chr8', 100, 150, False, False, False),
                samplot.PairedEnd('chr8', 170, 180, True, False, False)]

        max_event, steps = samplot.get_pairs_plan(ranges, pairs[hp])

        self.assertEqual(80, max_event)
        self.assertEqual(2, len(steps))
    #}}}
#}}}

#{{{ class Test_linked(unittest.TestCase):
class Test_linked(unittest.TestCase):
    #{{{def test_get_split_insert_size(self):
    def test_get_linked_plan(self):

        ranges = [samplot.genome_interval('chr8', 100, 200),
                  samplot.genome_interval('chr8', 600, 800) ]

        pairs = {}
        hp = 0
        pairs[hp] = {}

        pairs[hp]['PE_1'] = [\
                samplot.PairedEnd('chr8', 100, 150, False, False, False),
                samplot.PairedEnd('chr8', 170, 180, True, False, False)]

        pairs[hp]['PE_2'] = [\
                samplot.PairedEnd('chr8', 110, 160, False, False, False),
                samplot.PairedEnd('chr8', 680, 690, True, False, False)]

        splits = {}
        splits[hp] = {}

        splits[hp]['SR_1'] = [\
                samplot.SplitRead('chr8', 155, 160, True, 0, False, False),
                samplot.SplitRead('chr8', 670, 675, True, 50, False, False)]

        linked_reads = {}
        linked_reads[hp] = {}

        MI = 5
        linked_reads[hp][MI] = [[],[]]

        linked_reads[hp][MI][0].append('PE_1')
        linked_reads[hp][MI][0].append('PE_2')
        linked_reads[hp][MI][1].append('SR_1')


        max_event, steps = samplot.get_linked_plan(ranges,
                                                   pairs[hp],
                                                   splits[hp],
                                                   linked_reads[hp],
                                                   MI)
        self.assertEqual(580, max_event)
        self.assertEqual(2, len(steps))
        self.assertEqual(2, len(steps[0].info['PAIR_STEPS']))
        self.assertEqual(1, len(steps[0].info['SPLIT_STEPS']))

        self.assertEqual(100, steps[0].start_pos.start)
        self.assertEqual(100, steps[0].start_pos.end)
        self.assertEqual(ranges[0].end, steps[0].end_pos.start)
        self.assertEqual(ranges[0].end, steps[0].end_pos.end)

        self.assertEqual(ranges[1].start, steps[1].start_pos.start)
        self.assertEqual(ranges[1].start, steps[1].start_pos.end)
        self.assertEqual(690, steps[1].end_pos.start)
        self.assertEqual(690, steps[1].end_pos.end)


        self.assertEqual(100,steps[0].info['PAIR_STEPS'][0].start_pos.start)
        self.assertEqual(100,steps[0].info['PAIR_STEPS'][0].start_pos.end)
        self.assertEqual(180,steps[0].info['PAIR_STEPS'][0].end_pos.start)
        self.assertEqual(180,steps[0].info['PAIR_STEPS'][0].end_pos.end)

        self.assertEqual(110,steps[0].info['PAIR_STEPS'][1].start_pos.start)
        self.assertEqual(110,steps[0].info['PAIR_STEPS'][1].start_pos.end)
        self.assertEqual(690,steps[0].info['PAIR_STEPS'][1].end_pos.start)
        self.assertEqual(690,steps[0].info['PAIR_STEPS'][1].end_pos.end)


        self.assertEqual(160,steps[0].info['SPLIT_STEPS'][0].start_pos.start)
        self.assertEqual(160,steps[0].info['SPLIT_STEPS'][0].start_pos.end)
        self.assertEqual(670,steps[0].info['SPLIT_STEPS'][0].end_pos.start)
        self.assertEqual(670,steps[0].info['SPLIT_STEPS'][0].end_pos.end)

    #}}}
#}}}

if __name__ == '__main__':
    unittest.main()
