library(yaml)

load_config <- function(config_file = "config/config.yaml") {
  cfg <- yaml::read_yaml(config_file)
  return(cfg)
}

# Example usage inside your scripts:
# cfg <- load_config()
# patient_id <- cfg$patient_id
# genome_fasta <- cfg$genome_fasta
