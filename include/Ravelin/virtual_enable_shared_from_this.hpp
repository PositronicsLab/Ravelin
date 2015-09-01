#ifndef _VIRTUAL_SHARED_FROM_THIS_HPP_
#define _VIRTUAL_SHARED_FROM_THIS_HPP_

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

namespace Ravelin {

struct virtual_enable_shared_from_this_base:
   boost::enable_shared_from_this<virtual_enable_shared_from_this_base> {
   virtual ~virtual_enable_shared_from_this_base() {}
};

template<typename T>
struct virtual_enable_shared_from_this:
virtual virtual_enable_shared_from_this_base {
   boost::shared_ptr<T> shared_from_this() {
      return boost::dynamic_pointer_cast<T>(
         virtual_enable_shared_from_this_base::shared_from_this());
   }

   boost::shared_ptr<const T> shared_from_this() const {
      return boost::dynamic_pointer_cast<const T>(
         virtual_enable_shared_from_this_base::shared_from_this());
   }
};

}
#endif

