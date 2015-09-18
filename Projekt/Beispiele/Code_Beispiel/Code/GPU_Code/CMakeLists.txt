# CMake Code-Verzeichnis

# Aufsetzen
cmake_minimum_required(VERSION 2.8)

# Preprocessor-Variable DIR_VECTOR
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -DDIR_CUDA=\"${CMAKE_CURRENT_SOURCE_DIR}\"")

# Erzeugen der Executables
cuda_add_executable(SMVP_DEVICE smvp_device.cu)
cuda_add_executable(SMVP_CPU smvp_cpu.cu)
cuda_add_executable(CSR_DEVICE csr_smvp_device.cu)
#cuda_add_executable(CUSPARSE_DEVICE cusparse_smvp_device.cu)
# CUSPARSE Support muss noch gemacht werden
# zunaechst brauche Ich eine Liste aller zusaetzlichen Packages


