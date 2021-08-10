#pragma once

//  Use traditional AAD of chapter 10 (false)
//      or expression templated (AADET) of chapter 15 (true)
#ifndef AADET
#define AADET true
#endif

// makes it possible to have multiple derivatives at once
#ifndef AADMULTIOUT
#define AADMULTIOUT false
#endif