"""
A Python 3 wrapper for the Genome Analysis Toolkit.
"""

class Gatk(object):
    """
    NAME
        gatk -- A Python class wrapper for the Genome Analysis Toolkit

    SYNOPSIS
        object = gatk()

    DESCRIPTION
        A class that wraps the methods and attributes of the Gneome Analysis
        toolkit.
    """

    def __init__(self,
                 gatk="GenomeAnalysisTK.jar",
                 reference="ref.fa",
                 java='/usr/bin/java',
                 tmpdir='./tmp',
                 logginglevel='INFO',
                 output_dir='.'):
        """
        This is the standard object initialization method that is automatically
        executed when an object is instantiated with this class.

        NOTE:
            Do not call this method directly.

        USAGE:
            object = gatk(gatk='/usr/bin/GenomeAnalysisTk.jar', reference='hg19.fa')

        INPUT:
             * gatk: full path to the GATK JAR file (required)
             * reference: full path to the reference genome used (required)
             * java: full path to the java program (default: /usr/bin/java)
             * tmpdir: full path to the temporary directory (default: ./tmp)
             * logginglevel: the level of logging to be used (default: INFO)
             * output_dir: the output directory to output files to (default: .)
        """
        self.java = java
        self.gatk = gatk
        self.reference = reference
        self.tmpdir = tmpdir
        self.logginglevel = logginglevel
        self.output_dir = output_dir

        return None

    def file_cleanup(self, files):
        """
        A method for cleaning up intermediate files.

        USAGE:
            object.file_cleanup(
                files
                )

        INPUT:
            * files: an array of files to cleanup

        OUTPUT:
            Returns a list containing a command for file cleanup.
        """
        program = 'rm -f'
        files_string = ' '.join([files])
        cmd = ' '.join([program, files_string])
        return {"command":cmd}


    def realigner_target_creator(self,
                                 bams,
                                 output,
                                 memory=20,
                                 threads=8,
                                 knownsites=None):
        """
        Method wrapper for GATK's RealignerTargetCreator subprogram.

        USAGE:
            object.realigner_target_creator(
                bams = [input1, input2, input3],
                output = <output realigner target file>,
                memory = <memory allocated to Java>,
                threads = <number of threads to use>,
                knownsites = [knownsite1, knownsite2, knownsite3],
                output_dir = <output directory>
                )

        INPUT:
            * bams:         array of input files (required)
            * output:       output realigner target filename (required)
            * memory:       amount of memory to allocate to the Java engine (default: 20)
            * threads:      number of threads to use (default: '8')
            * knownsites:   a list of known indel sites (default: [])
            * output_dir:   path to the output directory to write the output to (default: .)

        OUTPUT:
            Returns a dictionary containing the command and output filename
        """
        gatk_memory = ''.join(['-Xmx', str(memory), 'g'])
        gatk_tmp = ''.join(['-Djava.io.tmpdir=', self.tmpdir])
        program = ' '.join([
            self.java,
            gatk_memory,
            gatk_tmp,
            '-jar', str(self.gatk)
            ])

        output_filename = '/'.join([self.output_dir, output])

        # build the option for input BAM files
        input_option = " ".join(['-I ' + input_bam for input_bam in bams])

        options = ' '.join([
            '-T RealignerTargetCreator',
            input_option,
            '-o', output_filename,
            '-R', self.reference,
            '-nt', str(threads),
            '-l', self.logginglevel
            ])

        # build the option for known indel sites
        if knownsites:
            knownsites_option = " ".join(['--known ' + site for site in knownsites])
            options = " ".join([options, knownsites_option])

        cmd = ' '.join([program, options])
        return {"command": cmd, "output":output}


    def indel_realigner(self,
                        bams,
                        output,
                        intervals,
                        memory=20):
        """
        Python wrapper for handling GATK's IndelRealigner sub-program.

        USAGE:
            object.indel_realigner(
                bams = [input1, input2, input3],
                output = <output BAM file>,
                intervals = <output from realigner target creator>,
                memory = <memory allocated to Java>
                )

        INPUT:
            * bam:          array of input BAM files (required)
            * output:       indel realigned output BAM file (required)
            * intervals:    output from the realigner target creator (required)
            * memory:       amount of memory to allocate to the Java engine (default: 20)

        OUTPUT:
             Returns a dictionary containing the command and output filename
        """
        java_memory = ''.join(['-Xmx', str(memory), 'g'])
        gatk_tmp = ''.join(['-Djava.io.tmpdir=', self.tmpdir])
        program = ' '.join([
            self.java,
            java_memory,
            gatk_tmp,
            '-jar', str(self.gatk)
            ])

        input_option = " ".join(['-I ' + input_bam for input_bam in bams])

        options = ' '.join([
            '-T IndelRealigner',
            input_option,
            '-o', output,
            '-targetIntervals', intervals,
            '-R', self.reference
            ])
        cmd = ' '.join([program, options])
        return {'command':cmd, 'output':output}


    def base_quality_recalibration(self,
                                   bam, output,
                                   memory=20,
                                   knownsites=None):
        """
        A method that wraps GATK's base quality recalibration suprogram.

        USAGE:
            object.base_quality_recalibration():
                INPUT

        INPUT:
            * bam: BAM file to process
            * output: name of file to write output to
            * memory: memory (in GB) to allocate to the Java engine (default: 26)

        OUTPUT:
             list of return values
        """
        java_memory = ''.join(['-Xmx', str(memory), 'g'])
        java_tmp = ''.join(['-Djava.io.tmpdir=', self.tmpdir])
        program = ' '.join([
            self.java,
            java_memory,
            java_tmp,
            "-jar", str(self.gatk)
            ])
        options = ' '.join([
            '-T', 'BaseRecalibrator',
            '-I', bam,
            '-o', output,
            '-R', self.reference,
            # '-rf', 'BadCigar',
            '-U ALLOW_N_CIGAR_READS'
            # '-cov ReadGroupCovariate',
            # '-cov ContextCovariate',
            # '-cov CycleCovariate',
            # '-cov QualityScoreCovariate'
            ])
        if knownsites:
            knownsites_option = " ".join(['--knownSites ' + site for site in knownsites])
            options = " ".join([options, knownsites_option])
        cmd = ' '.join([program, options])
        return {'command':cmd, 'output':output}


    def print_reads(self, bam, output, reference, recaldata):
        """
        A method that wraps GATK's PrintReads subprogram.

        USAGE:
            object.print_reads(
                 arg1
                )

        INPUT:
            *  arg1:  description of arg1

        OUTPUT:
            A BAM file containing base quality recalibrated BAM files.
        """
        pass
        # java_memory = ''.join([memory, 'g'])
        # program = ' '.join([
        #     '-nct', threads,
        #     '-I', input,
        #     '-R', self.reference,
        #     '-BQSR', recaldata,
        #     '-o', recal,
        #     '-rf BadCigar',
        #     coverage
        #     ])

        #return {"cmd": cmd, "output": output}
