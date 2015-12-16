ff = runtests('std_tests')
assert(all([ff(:).Passed]))
display('All tests passed!')
