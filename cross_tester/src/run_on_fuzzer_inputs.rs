fn read_inputs(dir: &str, ext: &str) -> Vec<Vec<u8>> {
    use std::io::Read;
    use std::fs::{self};
    use std::path::Path;
    use std::fs::File;

    let dir = Path::new(dir);
    assert!(dir.is_dir());
    let mut results = vec![];
    for entry in fs::read_dir(dir).expect("must read the directory") {
        let entry = entry.expect("directory should contain files");
        let path = entry.path();
        if path.is_dir() {
            continue
        } else {
            let extension = path.extension();
            if extension.is_none() {
                if ext != "" {
                    continue
                }
            } else {
                let extension = extension.unwrap();
                if extension != ext {
                    continue
                }
            }
        }
        let mut buffer = Vec::new();
        let mut f = File::open(path).expect("must open file");
        f.read_to_end(&mut buffer).expect("must read bytes from file");
        results.push(buffer);
    }
    
    results
}

#[test]
fn cross_check_on_honggfuzz() {
    use super::run;

    let path = "../honggfuzz/hfuzz_workspace/fuzz_target_compare/";
    let ext = "fuzz";
    let inputs = read_inputs(path, ext);
    println!("Running on {} crash inputs", inputs.len());
    for input in inputs.iter() {
        run(&input[..]);
    }
}

#[test]
fn cross_check_on_libfuzzer() {
    use super::run;

    let path = "../fuzz/artifacts/fuzz_target_compare/";
    let ext = "";
    let inputs = read_inputs(path, ext);
    println!("Running on {} crash inputs", inputs.len());
    for input in inputs.iter() {
        run(&input[..]);
    }
}