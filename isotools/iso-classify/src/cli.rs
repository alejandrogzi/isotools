use clap::{Parser, Subcommand};
use config::ArgCheck;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    #[command(subcommand)]
    pub command: SubArgs,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,
}

impl Args {}

#[derive(Debug, Subcommand)]
pub enum SubArgs {
    #[command(name = "intron")]
    Intron {
        #[command(flatten)]
        args: IntronArgs,
    },
    #[command(name = "exon")]
    Exon {
        #[command(flatten)]
        args: ExonArgs,
    },
}

#[derive(Debug, Parser)]
pub struct IntronArgs {
    #[arg(
        short = 'i',
        long = "iso",
        required = true,
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Paths to IsoSeq's BED12 file(s) delimited by comma"
    )]
    pub iso: Vec<PathBuf>,

    #[arg(
        short = 'b',
        long = "blacklist",
        required = false,
        value_name = "PATH",
        value_delimiter = ',',
        num_args = 1..,
        help = "Path to BED4 file with blacklisted introns"
    )]
    pub blacklist: Vec<PathBuf>,

    #[arg(
        short = 'w',
        long = "bigwig",
        required = false,
        value_name = "PATH",
        num_args = 1,
        help = "Path to spliceAI directory [will asume 2 files per strand: acceptor and donor .bw]"
    )]
    pub spliceai: Option<PathBuf>,

    #[arg(
        long = "twobit",
        required = false,
        value_name = "PATH",
        num_args = 1..,
        help = "Path to genome 2bit file"
    )]
    pub twobit: Option<PathBuf>,

    #[arg(
        short = 't',
        long = "toga",
        required = false,
        value_name = "PATH",
        num_args = 1..,
        help = "Path to TOGA annotation .bed file"
    )]
    pub toga: Option<PathBuf>,
}

impl ArgCheck for IntronArgs {
    fn get_blacklist(&self) -> &Vec<PathBuf> {
        &self.blacklist
    }

    fn get_ref(&self) -> &Vec<PathBuf> {
        &self.iso
    }

    fn get_query(&self) -> &Vec<PathBuf> {
        todo!()
    }
}

#[derive(Debug, Parser)]
pub struct ExonArgs {}
