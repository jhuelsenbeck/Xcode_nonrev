#include <iostream>
#include "Alignment.h"
#include "Mcmc.h"
#include "Model.h"
#include "RandomVariable.h"
#include "Settings.h"

#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 1300)
#   include <immintrin.h>
#else
#   include <stdint.h>
#   if defined(_MSC_VER)
#       include <intrin.h>
#   endif
int check_xcr0_zmm(void);
void run_cpuid(uint32_t eax, uint32_t ecx, uint32_t* abcd);
#endif

void checkConfiguration(void);
int hasIntelKnlFeatures(void);
void printHeader(void);

/*
 To do: 1. We might want to speed up the calculations further still by compressing the data.
        2. We need to add stepping stone integration so we can calculate marginal likelihoods.
        3. We may want to do analyses on a fixed tree, so we should add a tree constructor
           which takes as input a Newick string.
        4. We should have a more sensible prior on branch lengths.
*/


int main(int argc, char* argv[]) {

    // print nifty header
    printHeader();
    
    // check computer configuration
    checkConfiguration();
    
    // get the user settings
    Settings mySettings(argc, argv);
    
    // make the random number object
    RandomVariable rv;

    // read the data matrix
    Alignment myAlignment( mySettings.getInputFileName() );
    //myAlignment.print();

    // set up the phylogenetic model
    Model myModel(&mySettings, &myAlignment, &rv);
    
    // run the Markov chain Monte Carlo analysis
    Mcmc myMcmc(&mySettings, &myModel, &rv, &myAlignment);
    myMcmc.runPowerPosterior();
    
    return 0;
}

void checkConfiguration(void) {

    if ( hasIntelKnlFeatures() )
        std::cout << "   * Computer has Intel Knights Landing features" << std::endl;
}

#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 1300)
#else
int check_xcr0_zmm(void) {

    uint32_t xcr0;
    uint32_t zmm_ymm_xmm = (7 << 5) | (1 << 2) | (1 << 1);
#   if defined(_MSC_VER)
    xcr0 = (uint32_t)_xgetbv(0);  /* min VS2010 SP1 compiler is required */
#   else
    __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx" );
#   endif
    return ((xcr0 & zmm_ymm_xmm) == zmm_ymm_xmm); /* check if xmm, zmm and zmm state are enabled in XCR0 */
}

void run_cpuid(uint32_t eax, uint32_t ecx, uint32_t* abcd) {

#   if defined(_MSC_VER)
    __cpuidex(abcd, eax, ecx);
#   else
    uint32_t ebx, edx;
#   if defined( __i386__ ) && defined ( __PIC__ )
    /* in case of PIC under 32-bit EBX cannot be clobbered */
    __asm__ ( "movl %%ebx, %%edi \n\t cpuid \n\t xchgl %%ebx, %%edi" : "=D" (ebx),
#   else
    __asm__ ( "cpuid" : "+b" (ebx),
#   endif
    "+a" (eax), "+c" (ecx), "=d" (edx) );
    abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;
#   endif
}
#endif

int hasIntelKnlFeatures(void) {

#   if defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 1300)

    const unsigned long knl_features = (_FEATURE_AVX512F | _FEATURE_AVX512ER | _FEATURE_AVX512PF | _FEATURE_AVX512CD );
    return _may_i_use_cpu_feature( knl_features );
    
#   else

    uint32_t abcd[4];
    uint32_t osxsave_mask = (1 << 27); // OSX.
    uint32_t avx2_bmi12_mask = (1 << 16) | // AVX-512F
                               (1 << 26) | // AVX-512PF
                               (1 << 27) | // AVX-512ER
                               (1 << 28);  // AVX-512CD
    run_cpuid( 1, 0, abcd );
    // step 1 - must ensure OS supports extended processor state management
    if ( (abcd[2] & osxsave_mask) != osxsave_mask )
        return 0;
    // step 2 - must ensure OS supports ZMM registers (and YMM, and XMM)
    if ( ! check_xcr0_zmm() )
        return 0;
    return 1;
    
#   endif
}

void printHeader(void) {

    std::cout << std::endl;
    std::cout << "   Kasey's First Phylogeny Program, version 1.0" << std::endl;
    std::cout << std::endl;
    std::cout << "   * John Huelsenbeck and Kasey Arzumanova" << std::endl;
    std::cout << "   * University of California, Berkeley" << std::endl;
    std::cout << "   * " << std::endl;
}
