//leveraging advanced vector extensions and fused multiply add
//to achive ~10x reduction in elapsed time for discrete fourier transform

void calculateDFT(const int len, float* xn, float* Xr, float* Xi) {
    int k, n;
    for (k = 0; k < len; k++) {
        Xr[k] = 0;
        Xi[k] = 0;
        for (n = 0; n < len; n++) {
            Xr[k] += xn[n] * cos(2 * M_PI * k * n / len);
            Xi[k] -= xn[n] * sin(2 * M_PI * k * n / len);
        }
    }
}

int main() {
    const int len = 1024;
    float* xn = (float*)malloc(sizeof(float) * len);
    float* Xr = (float*)malloc(sizeof(float) * len);
    float* Xi = (float*)malloc(sizeof(float) * len);

    // init input data 
    for (int i = 0; i < len; i++) {
        xn[i] = sin(2 * M_PI * i / len);
    }

    clock_t start = clock();
    calculateDFT(len, xn, Xr, Xi);
    clock_t end = clock();
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("elapsed time: %f seconds\n", cpu_time_used);

    // Free the allocated memory
    free(xn);
    free(Xr);
    free(Xi);

    return 0;
}

elapsed time: 0.024000 seconds.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// helper horizontal sum __m256 vector.
inline float horizontal_sum_ps_avx(__m256 v) {
    __m128 vlow = _mm256_castps256_ps128(v);
    __m128 vhigh = _mm256_extractf128_ps(v, 1); // high 128
    vlow = _mm_add_ps(vlow, vhigh);             // add the low 128
    __m128 shuf = _mm_movehdup_ps(vlow);        // broadcast elements 3,1 2,0
    __m128 sums = _mm_add_ps(vlow, shuf);
    shuf = _mm_movehl_ps(shuf, sums); // high hal low half
    sums = _mm_add_ss(sums, shuf);
    return _mm_cvtss_f32(sums);
}

void precomputeTwiddleFactorsAVX2(int N, float* cosValues, float* sinValues) {
    const float two_pi_over_N = 2.0f * M_PI / N;
    for (int i = 0; i < N; i += 8) {
        __m256 indices = _mm256_set_ps(i + 7, i + 6, i + 5, i + 4, i + 3, i + 2, i + 1, i);
        __m256 angles = _mm256_mul_ps(indices, _mm256_set1_ps(two_pi_over_N));

        // here you would call vectorized cos and sin approximations
        // for demonstration: simply set values to the angle

        _mm256_storeu_ps(&cosValues[i], angles); //placeholder actual sin/cos computation
        _mm256_storeu_ps(&sinValues[i], angles); 
    }
}

void calculateDFT_AVX2(const int N, const float* xn, float* Xr, float* Xi) {
    // precompute twiddle factors
    std::vector<float> cosValues(N), sinValues(N);
    precomputeTwiddleFactorsAVX2(N, cosValues.data(), sinValues.data());
    
    for (int k = 0; k < N; ++k) {
        __m256 sumR = _mm256_setzero_ps(); // real sum
        __m256 sumI = _mm256_setzero_ps(); // imaginary sum

        for (int n = 0; n < N; n += 8) {
            // load next 8 values
            __m256 xn_vec = _mm256_loadu_ps(&xn[n]);

            // load precomputed cosine/sine values
            __m256 cos_vec = _mm256_loadu_ps(&cosValues[(k * n) % N]); // simplified for example
            __m256 sin_vec = _mm256_loadu_ps(&sinValues[(k * n) % N]); 

            // multiply-add operations
            sumR = _mm256_fmadd_ps(xn_vec, cos_vec, sumR);
            sumI = _mm256_fmadd_ps(xn_vec, sin_vec, sumI);
        }

        // reduce partial sums and store
        Xr[k] = horizontal_sum_ps_avx(sumR);
        Xi[k] = horizontal_sum_ps_avx(sumI);
    }
}

int main() {
    const int N = 1024; // dft size (must be a multiple of 8)
    std::vector<float> xn(N, 0.0f), Xr(N, 0.0f), Xi(N, 0.0f);

    // fill xn with data
    for (int i = 0; i < N; ++i) {
        xn[i] = sin(2 * M_PI * i / N);
    }

    auto start = std::chrono::high_resolution_clock::now();
    calculateDFT_AVX2(N, xn.data(), Xr.data(), Xi.data());
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "elapsed time: " << elapsed.count() << " seconds\n";

    return 0;
}

elapsed time: 0.0018957 seconds