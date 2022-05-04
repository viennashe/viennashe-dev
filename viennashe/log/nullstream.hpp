#ifndef VIENNASHE_LOG_NULLSTREAM_HPP
#define VIENNASHE_LOG_NULLSTREAM_HPP

/* ============================================================================
   Copyright (c) 2011-2022, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
     ViennaSHE - The Vienna Spherical Harmonics Expansion Boltzmann Solver
                            -----------------

                    http://viennashe.sourceforge.net/

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */

#include <iosfwd>
#include <memory>

/** @file viennashe/log/nullstream.hpp
    @brief Provides a streamer which discards all inputs. Similar in effect to a redirection to /dev/null. Used if logging is disabled.
*/

namespace viennashe
{
  namespace log
  {

    /** @brief Streaming class which only provides operator<< discarding the right hand side */
    struct nullstream {};

    /**
     * @brief Singleton factory for nullstream
     * @return Returns a non-const reference to the nullstream
     */
    inline nullstream & get_nullstream()
    {
      // not exactly singleton, but enough for our purpose
      //static std::auto_ptr<nullstream> stnull(new nullstream());
      static std::unique_ptr<nullstream> stnull(new nullstream());
      return *stnull;
    }

    /**
     * @brief Generic shift left operator. Throws anything it gets away.
     * @return The same object as get_nullstream()
     */
    template <typename T>
    nullstream & operator<<(nullstream & s, T const &)           { return s; }

    /** @brief Specialization of the generic template function. Returns same as get_nullstream(). Does nothing. */
    inline nullstream & operator<<(nullstream & s, const std::string &) { return s; }
    /** @brief Specialization of the generic template function. Returns same as get_nullstream(). Does nothing. */
    inline nullstream & operator<<(nullstream & s, const char * )       { return s; }

    /** @brief Specialization of the generic template function for std::endl. Returns same as get_nullstream(). Does nothing. */
    inline nullstream & operator<<(nullstream & s, std::ostream &(std::ostream&))  { return s; }


  } // log
} // viennashe


#endif

