#Variable for lauch snakemake
cluster_config='config_cluster.yaml'
snakemake -s annotation_pipeline.snake --latency-wait 555555 --jobs 100 --cluster "{cluster.scheduler} {cluster.queue} {cluster.export_env} {cluster.cwd} {cluster.mem} {cluster.n_cpu} {threads} " --cluster-config ${cluster_config} 
