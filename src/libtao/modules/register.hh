#ifndef tao_modules_register_hh
#define tao_modules_register_hh

#include "libtao/base/factory.hh"
#include "modules.hh"

namespace tao {

   template< class Backend >
   void
   register_modules( tao::factory<Backend>& fact )
   {
      fact.register_module( "light-cone", tao::modules::lightcone<Backend>::factory );
      fact.register_module( "hdf5", modules::hdf5<Backend>::factory );
   }

}

#endif
