extern crate hex;

#[test]
fn run_on_input() {
    let filename = "slow-unit-73a3ff99982d484fed518f581ba518394962c29f";
    use std::time::Instant;
    use std::io::Read;
    use std::fs::File;
    let mut file = File::open(&format!("fuzz/artifacts/fuzz_target_api/{}", filename)).expect("must open");
    let mut input_data = vec![];
    file.read_to_end(&mut input_data).expect("must read");
    assert!(input_data.len() != 0);
    let now = Instant::now();
    let result = crate::public_interface::API::run(&input_data[..]);
    let elapsed = now.elapsed().as_micros();
    let gas_estimate = crate::gas_meter::GasMeter::meter(&input_data[..]);
    if result.is_err() {
        println!("Api call failed in {} micros", elapsed);
        if gas_estimate.is_ok() {
            println!("Gas estimate was {}", gas_estimate.unwrap());
        } else {
            println!("Gas estimate failed");
        }
    } else {
        println!("Api call was ok in {} micros", elapsed);
        if gas_estimate.is_ok() {
            println!("Gas estimate was {}", gas_estimate.unwrap());
        } else {
            println!("Gas estimate failed");
        }
    }
}