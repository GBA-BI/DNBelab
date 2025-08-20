version 1.0

workflow DNBelab_C_Workflow {
    input {
        String sample_name
        Array[File] CDNAfastq1
        Array[File] CDNAfastq2
        Array[File] oligofastq1
        Array[File] oligofastq2
        Int cpu
        String memory = "128 GB"
        String disk = "1000 GB"
        String genome_name = "GRCm38"
    }

    String ref_index = if genome_name == "GRCm38" then "s3://gznl-open-data-references/BGI_SC/GRCm38.tar.gz"
                      else "s3://gznl-open-data-references/BGI_SC/GRCh38.tar.gz"

    call DNBelabanalysis {
        input:
            sample_name = sample_name,
            CDNAfastq1 = CDNAfastq1,
            CDNAfastq2 = CDNAfastq2,
            oligofastq1 = oligofastq1,
            oligofastq2 = oligofastq2,
            genome_name = genome_name,
            cpu = cpu,
            memory = memory,
            disk = disk,
            ref_index = ref_index
    }
    
    output {
        File report_html = DNBelabanalysis.report_html
        File anno_sorted_bam = DNBelabanalysis.anno_sorted_bam
        File anno_sorted_bam_bai = DNBelabanalysis.anno_sorted_bam_bai
        File filter_feature_h5ad = DNBelabanalysis.filter_feature_h5ad
        File metrics_summary = DNBelabanalysis.metrics_summary
        File singlecell_file = DNBelabanalysis.singlecell_file
        File filter_matrix = DNBelabanalysis.filter_matrix
        File raw_matrix = DNBelabanalysis.raw_matrix
        File RNAvelocity_matrix = DNBelabanalysis.RNAvelocity_matrix
        File splice_matrix = DNBelabanalysis.splice_matrix
        File all_results = DNBelabanalysis.all_results
        File output_dir = DNBelabanalysis.output_dir
    }
}

task DNBelabanalysis {
    input {
        String sample_name
        Array[File] CDNAfastq1
        Array[File] CDNAfastq2
        Array[File] oligofastq1
        Array[File] oligofastq2
        File ref_index
        Int cpu
        String memory
        String disk
        String genome_name
    }
    
    command <<<
        pigz -p 24 -dc ~{ref_index} | tar -xvf - -C ${PWD}
        mkdir /opt/~{genome_name}
        mv ~{genome_name}/* /opt/~{genome_name}
        /opt/dnbc4tools2.1.3/dnbc4tools rna run --cDNAfastq1 ~{sep="," CDNAfastq1} --cDNAfastq2 ~{sep="," CDNAfastq2} --oligofastq1 ~{sep="," oligofastq1} --oligofastq2 ~{sep="," oligofastq2} --genomeDir /opt/~{genome_name} --name ~{sample_name} --threads 24
        # 先打包原始输出目录
        tar -zcvf ~{sample_name}.tar.gz ~{sample_name}
        
        # 打包output目录
        tar -zcvf ~{sample_name}_output.tar.gz ~{sample_name}/output
        
        mv ./log ./log_old
        mv ~{sample_name}/output/* ./
        tar -zcvf raw_matrix.tar.gz raw_matrix
        tar -zcvf filter_matrix.tar.gz filter_matrix
        mv attachment/* ./
        tar -zcvf RNAvelocity_matrix.tar.gz RNAvelocity_matrix
        tar -zcvf splice_matrix.tar.gz splice_matrix
    >>>
    
    runtime {
        docker: "registry-vpc.miracle.ac.cn/gznl/dnbctools:v1"
        cpu: cpu
        memory: memory
        disk: disk
    }
    
    output {
        File report_html = "~{sample_name}_scRNA_report.html"
        File anno_sorted_bam = "anno_decon_sorted.bam"
        File anno_sorted_bam_bai = "anno_decon_sorted.bam.bai"
        File filter_feature_h5ad = "filter_feature.h5ad"
        File metrics_summary = "metrics_summary.xls"
        File singlecell_file = "singlecell.csv"
        File filter_matrix = "filter_matrix.tar.gz"
        File raw_matrix = "raw_matrix.tar.gz"
        File RNAvelocity_matrix = "RNAvelocity_matrix.tar.gz"
        File splice_matrix = "splice_matrix.tar.gz"
        File all_results = "~{sample_name}.tar.gz"
        File output_dir = "~{sample_name}_output.tar.gz"
    }
}
