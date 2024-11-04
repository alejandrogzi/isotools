use hashbrown::HashSet;
use packbed::par_reader;
use rayon::iter::ParallelIterator;
use rayon::str::ParallelString;
use std::path::PathBuf;

pub fn unpack_blacklist<'a>(paths: Vec<PathBuf>) -> Option<HashSet<String>> {
    if paths.is_empty() {
        return None;
    }

    let contents = par_reader(paths).unwrap();
    let tracks = contents
        .par_lines()
        .filter_map(|line| {
            if line.is_empty() {
                return None;
            };

            Some(line.to_string())
        })
        .collect::<HashSet<String>>();

    Some(tracks)
}
