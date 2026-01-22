/*
 * Test program to call the c_perform_calculation function from the multiply_numbers.dylib library.
 *
 * Instructions:
 * 1. Ensure the multiply_numbers.dylib library is built in build/lib/ (run julia multiply_numbers.jl).
 * 2. Compile this program: gcc test_c_call_multiply_numbers.c -o test_c_call_multiply_numbers
 * 3. Run the program: ./test_c_call_multiply_numbers
 *    Expected output: Result: 450
 *
 * This demonstrates calling a Julia @ccallable function from C using dlopen/dlsym.
 */

#include <stdio.h>
#include <dlfcn.h>

int main() {
    void* lib = dlopen("build/lib/multiply_numbers.dylib", RTLD_LAZY);
    if (!lib) {
        fprintf(stderr, "Error loading library: %s\n", dlerror());
        return 1;
    }
    
    long long (*c_perform_calculation)(long long*, int) = dlsym(lib, "c_perform_calculation");
    if (!c_perform_calculation) {
        fprintf(stderr, "Error finding symbol: %s\n", dlerror());
        dlclose(lib);
        return 1;
    }
    
    long long numbers[] = {5, 9, 10};
    long long result = c_perform_calculation(numbers, 3);
    printf("Result: %lld\n", result);  // Should print 450
    
    dlclose(lib);
    return 0;
}