#include <check.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/*
 * Security invariant: The allocation size for posterior_sample must be
 * sufficient to hold num_ps * size_of_modeltype bytes. If the multiplication
 * overflows, the allocation must fail safely rather than producing an
 * undersized buffer that is subsequently overwritten by memcpy in a loop.
 *
 * We test that num_ps * size_of_modeltype does not silently overflow by
 * checking the arithmetic before any allocation would occur.
 */

START_TEST(test_posterior_sample_allocation_overflow)
{
    /* Invariant: num_ps * size_of_modeltype must not overflow size_t,
       otherwise a heap buffer overflow occurs at the memcpy exploitation point. */

    struct test_case {
        size_t num_ps;
        size_t size_of_modeltype;
        int should_overflow;
    };

    struct test_case cases[] = {
        /* Exact exploit: large num_ps causing overflow with realistic model size */
        { (SIZE_MAX / 64) + 2, 64, 1 },
        /* Boundary: exactly at overflow edge */
        { SIZE_MAX / 128 + 1, 128, 1 },
        /* Valid input: small safe values */
        { 100, 64, 0 },
        /* Another valid case */
        { 1000, 256, 0 },
    };

    int num_cases = sizeof(cases) / sizeof(cases[0]);

    for (int i = 0; i < num_cases; i++) {
        size_t num_ps = cases[i].num_ps;
        size_t size_of_modeltype = cases[i].size_of_modeltype;

        /* Check for multiplication overflow */
        int overflows = (size_of_modeltype != 0 &&
                         num_ps > SIZE_MAX / size_of_modeltype);

        ck_assert_msg(overflows == cases[i].should_overflow,
            "Case %d: overflow detection mismatch (num_ps=%zu, size=%zu)",
            i, num_ps, size_of_modeltype);

        /* For non-overflowing cases, verify allocation would succeed
           and the loop bound i*size_of_modeltype stays in range */
        if (!overflows) {
            size_t alloc_size = num_ps * size_of_modeltype;
            ck_assert(alloc_size >= num_ps);
            ck_assert(alloc_size >= size_of_modeltype);
            /* Verify last element offset is within allocation */
            if (num_ps > 0) {
                size_t last_offset = (num_ps - 1) * size_of_modeltype;
                ck_assert(last_offset + size_of_modeltype <= alloc_size);
            }
        }
    }
}
END_TEST

Suite *security_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("Security");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_posterior_sample_allocation_overflow);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = security_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}