#ifndef LEDGER_HH
#define LEDGER_HH

#include <cstddef>

void ledger_init();
void ledger_alloc(const void* host, std::size_t size);
void ledger_free(const void* host);
void ledger_copyToDevice(const void* host);
void ledger_copyToHost(void* host);
void* ledger_lookup(const void* host);

#endif
