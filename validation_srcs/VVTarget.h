// include VVTARGET according to defined macro
//

#if defined(VVTARGET_TestOklahoma)
#include "TestOklahoma.h"
using VVTarget = TestOklahoma;

#elif defined(VVTARGET_TestCavityFlow)
#include "TestCavityFlow.h"
using VVTarget = TestCavityFlow;

#elif defined(VVTARGET_TestFlowAroundCube)
#include "TestFlowAroundCube.h"
using VVTarget = TestFlowAroundCube;

#else
#error no VVTARGET_* is defined

#endif
