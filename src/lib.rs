use pyo3::prelude::{Python, PyModule, PyErr};

#[no_mangle]
pub extern "C" fn multiply(factor0: f32, factor1: f32) -> f32 {
    let product = Python::with_gil(|py| -> Result<f32, PyErr> {
        println!("Message from Rust.");
        let python_module = PyModule::import(py, "eesunhong.python_module")?;
        let args = (factor0, factor1);
        let product: f32 = python_module.getattr("multiply")?.call1(args,)?.extract().unwrap();
        println!("Product in Rust is {}.", product);
        Ok(product)
    }).unwrap();
    product
}
