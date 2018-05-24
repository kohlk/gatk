# CRAM to filtered VCF WDL
workflow Cram2TrainedModel {
    File input_cram
    File reference_fasta
    File reference_dict
    File reference_fasta_index
    File truth_vcf
    File truth_vcf_index
    File truth_bed
    String output_prefix
    String tensor_type
    Int epochs
    File calling_intervals
    Int scatter_count
    String extra_args

    # Runtime parameters
    File? gatk_override
    String gatk_docker
    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu

    call CramToBam {
        input:
          reference_fasta = reference_fasta,
          reference_dict = reference_dict,
          reference_fasta_index = reference_fasta_index,
          cram_file = input_cram,
          output_prefix = output_prefix,
          disk_size = disk_space_gb
    }

    call SplitIntervals {
        input:
            gatk_override = gatk_override,
            scatter_count = scatter_count,
            intervals = calling_intervals,
            ref_fasta = reference_fasta,
            ref_dict = reference_dict,
            ref_fai = reference_fasta_index,
            preemptible_attempts = preemptible_attempts
    }

    scatter (calling_interval in SplitIntervals.interval_files) {

        call RunHC4 {
            input:
                input_bam = CramToBam.output_bam,
                input_bam_index = CramToBam.output_bam_index,
                reference_fasta = reference_fasta,
                reference_dict = reference_dict,
                reference_fasta_index = reference_fasta_index,
                output_prefix = output_prefix,
                interval_list = calling_interval,
                gatk_docker = gatk_docker,
                gatk_override = gatk_override,
                preemptible_attempts = preemptible_attempts,
                extra_args = extra_args,
                disk_space_gb = disk_space_gb
        }

        call WriteTensors {
            input:
                input_vcf = RunHC4.raw_vcf,
                input_vcf_index = RunHC4.raw_vcf_index,
                input_bam = RunHC4.bamout,
                input_bam_index = RunHC4.bamout_index,
                truth_vcf = truth_vcf,
                truth_vcf_index = truth_vcf_index,
                truth_bed = truth_bed,
                reference_fasta = reference_fasta,
                reference_dict = reference_dict,
                reference_fasta_index = reference_fasta_index,
                output_prefix = output_prefix,
                interval_list = calling_interval,
                gatk_jar = gatk_override,
                disk_size = disk_space_gb,
                tensor_type = tensor_type
        }
    }

    call MergeVCFs as MergeVCF_CNN {
        input:
            input_vcfs = RunHC4.raw_vcf,
            output_vcf_name = output_prefix,
            gatk_override = gatk_override,
            gatk_docker = gatk_docker
    }

    call SamtoolsMergeBAMs {
        input:
            input_bams = RunHC4.bamout,
            output_prefix = output_prefix
    }

    call TrainModel {
        input:
            tar_tensors = WriteTensors.tensors,
            output_prefix = output_prefix,
            tensor_type = tensor_type,
            gatk_jar = gatk_override,
            disk_size = disk_space_gb,
            epochs = epochs
    }

    output {
        MergeVCF_HC4.*
        SamtoolsMergeBAMs.*
        TrainModel.*
    }

}

task CramToBam {
  File reference_fasta
  File reference_fasta_index
  File reference_dict
  File cram_file
  String output_prefix

  Int disk_size

command <<<
  ls -ltr ${cram_file} ${reference_fasta} &&
  echo "ls (1): complete" &&
  samtools view -h -T ${reference_fasta} ${cram_file} |
  samtools view -b -o ${output_prefix}.bam - &&
  echo "samtools view: complete" &&
  ls -ltr . &&
  echo "ls (2): complete" &&
  samtools index -b ${output_prefix}.bam &&
  echo "samtools index: complete" &&
  ls -ltr . &&
  echo "ls (3): complete" &&
  mv ${output_prefix}.bam.bai ${output_prefix}.bai &&
  echo "mv: complete" &&
  ls -ltr . &&
  echo "ls (4): complete"
  >>>
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.1.1"
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "${output_prefix}.bam"
    File output_bam_index = "${output_prefix}.bai"
  }
}


task RunHC4 {
    File input_bam
    File input_bam_index
    File reference_fasta
    File reference_dict
    File reference_fasta_index
    String output_prefix
    File interval_list
    File gatk_jar
    String extra_args
    Int disk_size

    command { 
        java -Djava.io.tmpdir=tmp -jar ${gatk_jar} \
        HaplotypeCaller \
        -R ${reference_fasta} \
        -I ${input_bam} \
        -O ${output_prefix}_hc4.vcf.gz \
        -L ${interval_list} \
        -bamout ${output_prefix}_bamout.bam \
        ${extra_args}
    }

    output {
        File bamout = "${output_prefix}_bamout.bam"
        File bamout_index = "${output_prefix}_bamout.bai"
        File raw_vcf = "${output_prefix}_hc4.vcf.gz"
        File raw_vcf_index = "${output_prefix}_hc4.vcf.gz.tbi"
    }
    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.1.1"
        memory: "3 GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task WriteTensors {
    File input_bam
    File input_bam_index
    File input_vcf
    File input_vcf_index
    File reference_fasta
    File reference_dict
    File reference_fasta_index
    File truth_vcf
    File truth_vcf_index
    File truth_bed
    String output_prefix
    String tensor_type
    File interval_list
    File gatk_jar
    Int disk_size

    command {
        mkdir "./tensors/"

        java -Djava.io.tmpdir=tmp -jar ${gatk_jar} \
        CNNVariantWriteTensors \
        -R ${reference_fasta} \
        -V ${input_vcf} \
        -truth-vcf ${truth_vcf} \
        -truth-bed ${truth_bed} \
        -tensor-type ${tensor_type} \
        -output-tensor-dir "./tensors/" \
        -bam-file ${input_bam}
        
        tar -czf "tensors.tar.gz" "./tensors/"
    }

    output {
      File tensors = "tensors.tar.gz"
    }
    runtime {
        docker: "samfriedman/p3"
        memory: "3 GB"
        disks: "local-disk " + disk_size + " HDD"
    }

}


task TrainModel {
    Array[File] tar_tensors
    String output_prefix
    String tensor_type
    File gatk_jar
    Int disk_size
    Int epochs
    command {
        for tensors in ${sep=' ' tar_tensors}  ; do
            tar -xzf $tensors 
        done

        java -Djava.io.tmpdir=tmp -jar ${gatk_jar} \
        CNNVariantTrain \
        -input-tensor-dir "./tensors/" \
        -model-name ${output_prefix} \
        -image-dir "./" \
        -tensor-type ${tensor_type} \
        -epochs ${epochs}
    }

    output {
        File model_json = "${output_prefix}.json"
        File model_hd5 = "${output_prefix}.hd5"
        File roc_png = "per_class_roc_${output_prefix}.png"
        File training_png = "metric_history_${output_prefix}.png"
    }

    runtime {
        docker: "samfriedman/p3"
        memory: "16 GB"
        disks: "local-disk " + disk_size + " HDD"

#      docker: "samfriedman/gpu"
#      gpuType: "nvidia-tesla-k80" 
#      gpuCount: 1 
#      zones: ["us-central1-c"]
#      memory: "16 GB"
#      disks: "local-disk 400 HDD"
#      bootDiskSizeGb: "16" 
    }
}


task MergeVCFs {
    # inputs
    Array[File] input_vcfs
    String output_prefix

    File? gatk_override

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 1000

    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}
        #gatk --java-options "-Xmx${command_mem}m" MergeVcfs \
        java "-Xmx${command_mem}m" -jar ${gatk_override} MergeVcfs \
        -I ${sep=' -I ' input_vcfs} -O "${output_prefix}_cnn_scored.vcf.gz"
    }

    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + " SSD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }

    output {
        File merged_vcf = "${output_prefix}_hc4.vcf.gz"
        File merged_vcf_index = "${output_prefix}_hc4.vcf.gz.tbi"
    }
}


task SamtoolsMergeBAMs {
    Array[File] input_bams
    String output_prefix
    command {
        samtools merge ${output_prefix}_bamout.bam ${sep=' ' input_bams}
        samtools index ${output_prefix}_bamout.bam ${output_prefix}_bamout.bai
    }

    output {
        File bamout = "${output_prefix}_bamout.bam"
        File bamout_index = "${output_prefix}_bamout.bai"
    }

  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.1.1"
    memory: "16 GB"
    disks: "local-disk 400 HDD"
  }    
}


task SplitIntervals {
    # inputs
    File? intervals
    File ref_fasta
    File ref_fai
    File ref_dict
    Int scatter_count
    String? split_intervals_extra_args

    File? gatk_override

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        mkdir interval-files
        #gatk --java-options "-Xmx${command_mem}m" SplitIntervals \
        java "-Xmx${command_mem}m" -jar ${gatk_override} \
            SplitIntervals \
            -R ${ref_fasta} \
            ${"-L " + intervals} \
            -scatter ${scatter_count} \
            -O interval-files \
            ${split_intervals_extra_args}
        cp interval-files/*.intervals .
    }

    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + " SSD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }

    output {
        Array[File] interval_files = glob("*.intervals")
    }
}
