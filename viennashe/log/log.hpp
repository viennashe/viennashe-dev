#ifndef VIENNASHE_LOG_LOG_HPP
#define VIENNASHE_LOG_LOG_HPP

/* ============================================================================
   Copyright (c) 2011-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
     ViennaSHE - The Vienna Spherical Harmonics Expansion Boltzmann Solver
                            -----------------

                    http://viennashe.sourceforge.net/

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */

#include <iosfwd>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <stdio.h>

#include "viennashe/log/nullstream.hpp"

// Use this to enable fancy output for BASH
//#define VIENNASHE_LOG_ENABLE_FANCY_BASH

/** @file viennashe/log/log.hpp
    @brief A logging facility providing fine-grained control over logging in ViennaSHE
*/



namespace viennashe
{
  /** @brief Namespace holding the logging facility */
  namespace log
  {

    /** @brief Defines various log-levels. Note that log levels are inclusive, i.e. log levels with larger number include log levels with smaller number */
    enum log_levels
    {
      logERROR = 0,       /// logger outputs errors only
      logWARNING,         /// logger outputs errors and warnings
      logINFO,            /// logger outputs errors, warnings and status information
      logDEBUG,           /// logger outputs lots of debug infos as well

      logNEVER /// Nothing is logged. Internal note: must be the last in the list!
    };

    /** @brief Namespace for implementation details within the viennashe::log namespace. Typically not of interest for the library user. */
    namespace detail
    {
      static log_levels modify_log_level(log_levels new_level, bool setter = false)
      {
        static log_levels current_level = logDEBUG; // DEFAULT LOG LEVEL IS DEBUG
        if (setter) current_level = new_level;
        return current_level;
      }
    }

    /** @brief Sets the global log level */
    inline void set_log_level(log_levels new_level) { viennashe::log::detail::modify_log_level(new_level, true); }
    /** @brief Getter for the global log level */
    inline log_levels log_level() { return viennashe::log::detail::modify_log_level(viennashe::log::logDEBUG, false); }

    /**
     * @brief The Main logger class. Assembles output lines and writes them to std::cout upon destruction.
     * @tparam enabled If false the logger class does the same as the nullstream ... nothing
     */
    template < bool enabled >
    class logger
    {
      private:
        typedef std::ostringstream CollectorStreamType;
        CollectorStreamType local_out;

        logger & operator=(const logger &) { return *this; }

        log_levels _messageLevel;

      public:

        explicit logger() : _messageLevel(log_level()) {    }

        logger(log_levels level) : _messageLevel(level) {   }

        /** @brief CTOR to log componentwise. Adds [component_name] to the start of every log-line */
        logger(const std::string & component_name) : _messageLevel(log_level())
        {
          this->get() << "[" << component_name << "] ";
        }

        /** @brief CTOR to log componentwise on a certain log-level. Adds [component_name] to the start of every log-line */
        logger(log_levels level, const std::string & component_name) : _messageLevel(level)
        {
          this->get() << "[" << component_name << "] ";
        }

        logger(const logger & r) : _messageLevel(r._messageLevel)
        {
          if ( do_log() ) get() << r.get(); // _messageLevel must have been copied already!
        }

        /** @brief Destructor. Does actually write the log-message to the output-stream (normally std::cout) */
        ~logger()
        {
#ifndef VIENNASHE_LOG_DISABLE
          if (this->do_log())
          {
  #ifdef VIENNASHE_LOG_ENABLE_FANCY_BASH
            switch(_messageLevel)
            {
              case log::logERROR:
              #ifdef VIENNASHE_LOG_ATOMIC
                fprintf(stderr, "\033[1;31m %s \033[0m", get().str().c_str());
                fflush(stderr);
              #else
                std::cout << get().str();
              #endif
                break;
              case log::logWARNING:
              #ifdef VIENNASHE_LOG_ATOMIC
                fprintf(stdout, "\033[1;33m %s \033[0m", get().str().c_str());
                fflush(stdout);
              #else
                std::cout << get().str();
              #endif
                break;
              case log::logDEBUG:
              #ifdef VIENNASHE_LOG_ATOMIC
                fprintf(stdout, "\033[1;32m %s \033[0m", get().str().c_str());
                fflush(stdout);
              #else
                std::cout << get().str();
              #endif
                break;
              default:
              #ifdef VIENNASHE_LOG_ATOMIC
                fprintf(stdout, "%s", get().str().c_str());
                fflush(stdout);
              #else
                std::cout << get().str();
              #endif
                break;
            }
  #else
            #ifdef VIENNASHE_LOG_ATOMIC
              fprintf(stdout, "%s", get().str().c_str());
              fflush(stdout);
            #else
              std::cout << get().str();
            #endif
  #endif
          }
#endif
        }

        /** @brief Returns true if the log-level is smaller than the globally set one. */
        bool do_log() const { return (this->_messageLevel <= log_level()); }

        /** @brief Returns the collector stream to collect the output. */
        CollectorStreamType       & get()       {  return local_out;   }
        /** @brief Returns the collector stream to collect the output. */
        CollectorStreamType const & get() const {  return local_out;   }

        /** @brief Generic shift left operator to print stuff via the logger. */
        template <typename T>
        CollectorStreamType & operator<<(const T & x )
        {
          if( do_log() ) { get() << x; }
          return get();
        }

        typedef std::ostream& (*ostream_manipulator)(std::ostream&);
        CollectorStreamType & operator<<(ostream_manipulator pf)
        {
          if( do_log() ) { get() << pf ; }
          return get();
        }

        CollectorStreamType & operator<<(std::string & x )
        {
          if( do_log() ) { get() << x; }
          return get();
        }

        CollectorStreamType & operator<<(const std::string & x )
        {
          if( do_log() ) { get() << x; }
          return get();
        }

        CollectorStreamType & operator<<(const char * x )
        {
          if( do_log() && x != 0) { get() << x; }
          return get();
        }

        CollectorStreamType & operator<<(logger & r )
        {
          if( do_log() ) { get() << r.get(); }
          return get();
        }
    };

    /** @brief Template specialization of the logger for the case enabled=false ... does nothing. Ensures that no
     runtime penalty is present if the logger is disabled via the template parameter. */
    template < >
    class logger<false>
    {
      private:

        logger & operator=(const logger &) { return *this; }

        log_levels _messageLevel;

      public:

        explicit logger() : _messageLevel(log_level()) {    }

        logger(log_levels level) : _messageLevel(level) {   }

        logger(const std::string &) : _messageLevel(log_level()) {  }

        logger(log_levels level, const std::string &) : _messageLevel(level) {  }

        logger(const logger & r) : _messageLevel(r._messageLevel) { }

        ~logger() {  }

        bool do_log() const { return false; }

        nullstream & get() {  return get_nullstream();   }
        nullstream const& get() const {  return get_nullstream();   }

        template <typename T>
        nullstream & operator<<(const T &) { return get_nullstream(); }

        typedef std::ostream& (*ostream_manipulator)(std::ostream&);
        nullstream & operator<<(ostream_manipulator) { return get_nullstream(); }

        nullstream & operator<<(std::string &) { return get_nullstream(); }
        nullstream & operator<<(const std::string &) { return get_nullstream(); }
        nullstream & operator<<(const char *) { return get_nullstream(); }
        nullstream & operator<<(logger &) { return get_nullstream(); }
    };

    /** @brief Implementation details of the logging facility */
    namespace detail
    {

#ifndef VIENNASHE_LOG_DISABLE

    template < typename KeyTypeT >
    logger<KeyTypeT::enabled> vlogT()
    {
      return logger<KeyTypeT::enabled>();
    }

    template < log_levels level, typename KeyTypeT >
    logger<KeyTypeT::enabled> vlogTL()
    {
      return logger<KeyTypeT::enabled>(level);
    }

    template < log_levels level >
    logger<true> vlogL( )
    {
      return logger<true>(level); // using the given level
    }

#else
    template < typename KeyTypeT >
    nullstream & vlogT()
    {
      return get_nullstream();
    }

    template < log_levels level, typename KeyTypeT >
    nullstream & vlogTL()
    {
      return get_nullstream();
    }

    template < log_levels level >
    nullstream & vlogL( )
    {
      return get_nullstream();
    }
#endif

    } // namespace detail

#ifndef VIENNASHE_LOG_DISABLE

    /** @brief Used to log errors. The logging level is logERROR */
    inline logger<true>  error() { return detail::vlogL<logERROR>();   }
    /** @brief Used to log warnings. The logging level is logWARNING */
    inline logger<true>  warn()  { return detail::vlogL<logWARNING>(); }
    /** @brief Used to log warnings. The logging level is logWARNING */
    inline logger<true>  warning()  { return detail::vlogL<logWARNING>(); }
    /** @brief Used to log infos. The logging level is logINFO */
    inline logger<true>  info()  { return detail::vlogL<logINFO>();    }
    /** @brief Used to log debug output. The logging level is logDEBUG */
    inline logger<true>  debug() { return detail::vlogL<logDEBUG>();   }



    /** @brief Used to log errors for a certain component.
     *         If KeyTypeT::enabled is false no output will be generated. The logging level is logERROR */
    template < typename KeyTypeT >
    logger<KeyTypeT::enabled>  error() { return detail::vlogTL<logERROR, KeyTypeT>();   }

    /** @brief Used to log warnings for a certain component.
     *         If KeyTypeT::enabled is false no output will be generated. The logging level is logWARNING */
    template < typename KeyTypeT >
    logger<KeyTypeT::enabled>  warn()  { return detail::vlogTL<logWARNING, KeyTypeT>(); }

    /** @brief Used to log warnings for a certain component.
     *         If KeyTypeT::enabled is false no output will be generated. The logging level is logWARNING */
    template < typename KeyTypeT >
    logger<KeyTypeT::enabled>  warning()  { return detail::vlogTL<logWARNING, KeyTypeT>(); }

    /** @brief Used to log infos for a certain component.
     *         If KeyTypeT::enabled is false no output will be generated. The logging level is logINFO */
    template < typename KeyTypeT >
    logger<KeyTypeT::enabled>  info()  { return detail::vlogTL<logINFO, KeyTypeT>();    }

    /** @brief Used to log debug output for a certain component.
     *         If KeyTypeT::enabled is false no output will be generated. The logging level is logDEBUG */
    template < typename KeyTypeT >
    logger<KeyTypeT::enabled>  debug() { return detail::vlogTL<logDEBUG, KeyTypeT>();   }

#else

    inline nullstream & error() { return detail::vlogL<logERROR>();   }
    inline nullstream & warn()  { return detail::vlogL<logWARNING>(); }
    inline nullstream & warning()  { return detail::vlogL<logWARNING>(); }
    inline nullstream & info()  { return detail::vlogL<logINFO>();    }
    inline nullstream & debug() { return detail::vlogL<logDEBUG>();   }

    template < typename KeyTypeT >
    nullstream & error() { return detail::vlogTL<logERROR, KeyTypeT>();   }

    template < typename KeyTypeT >
    nullstream & warn()  { return detail::vlogTL<logWARNING, KeyTypeT>(); }

    template < typename KeyTypeT >
    nullstream & warning()  { return detail::vlogTL<logWARNING, KeyTypeT>(); }

    template < typename KeyTypeT >
    nullstream & info()  { return detail::vlogTL<logINFO, KeyTypeT>();    }

    template < typename KeyTypeT >
    nullstream & debug() { return detail::vlogTL<logDEBUG, KeyTypeT>();   }

#endif

  } // log
} // viennashe


#endif

