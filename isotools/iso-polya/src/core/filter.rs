use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};

use crate::{cli::FilterArgs, utils::get_assets_dir};

pub const PARA: &str = "para";
pub const ASSETS: &str = "assets";
pub const JOBLIST: &str = "joblist";
pub const FILTER_MINIMAP: &str = "filterMinimapQuality.perl";

pub fn filter_minimap(args: FilterArgs) -> Result<(), Box<dyn std::error::Error>> {
    let joblist = create_job_filter_minimap(args.clone());

    if args.para {
        let mut additional_args = vec![];
        if let Some(queue) = args.queue {
            additional_args.push(format!("-q {}", queue));
        }

        if let Some(mem) = args.mem {
            additional_args.push(format!("-memoryMb {}", mem));
        }

        let code = std::process::Command::new(PARA)
            .arg("make")
            .arg("filter_minimap")
            .arg(joblist)
            .args(additional_args)
            .output()
            .expect("ERROR: Failed to submit job");

        if !code.status.success() {
            let err = String::from_utf8_lossy(&code.stderr);
            log::error!("ERROR: Job failed to submit! {}", err);
            std::process::exit(1);
        } else {
            log::info!("SUCCESS: Jobs submitted!");
        }

        log::info!("INFO: Jobs finished successfully!");
    } else {
        log::info!(
            "{}",
            format!(
                "INFO: Only writing joblist to: {}",
                joblist.to_string_lossy()
            )
        );
    }

    Ok(())
}

fn create_job_filter_minimap(args: FilterArgs) -> PathBuf {
    let joblist = PathBuf::from(JOBLIST);
    let file = File::create(joblist.clone()).expect("Failed to create joblist");
    let mut writer = BufWriter::new(file);
    let executable = get_assets_dir().join(FILTER_MINIMAP);

    let mut cmd = format!(
        "perl {} {} -perID {} -clip3 {} -clip5 {} -P2P {} -emitA {}",
        executable.to_string_lossy(),
        args.sam
            .get(0)
            .expect("ERROR: No sam file provided")
            .to_string_lossy(),
        args.per_id,
        args.clip3,
        args.clip5,
        args.p2p,
        args.emit_a
    );

    if args.stat {
        cmd.push_str(" -statFile");
    }

    if args.keep {
        cmd.push_str(" -keepBad5Prime");
    }

    if let Some(suffix) = args.suffix {
        let suffix = format!(" -polyAReadSuffix {}", suffix);
        cmd.push_str(&suffix);
    }

    let _ = writer.write_all(cmd.as_bytes());
    let _ = writer.write_all(b"\n");

    let _ = writer.flush();

    return joblist;
}
