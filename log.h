#include <functional>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <utility>
#include <cstdio>
#include <ctime>
#include <chrono>
#include <thread>

#define FATAL      0
#define ERROR      1
#define WARN       2
#define INFO       3
#define DEBUG      4
#define TRACE      5

// Get current time, format is HH:mm:ss
const std::string currentTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[10];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%X", &tstruct);
    return buf;
}

typedef std::function<std::string()> prefixer;

// http://stackoverflow.com/questions/27336335/c-cout-with-prefix
class prefixbuf
    : public std::streambuf
{
protected:
  std::streambuf* sbuf;
private:
  prefixer m_prefix;
  bool     need_prefix;
public:
  int sync() {
    return this->sbuf->pubsync();
  }
  int overflow(int c) {
    if (c != std::char_traits<char>::eof()) {
      std::string prefix = m_prefix();
        if (this->need_prefix
            && !prefix.empty()
            && prefix.size() != this->sbuf->sputn(&prefix[0], prefix.size())) {
          return std::char_traits<char>::eof();
        }
      this->need_prefix = c == '\n';
    }
    return this->sbuf->sputc(c);
  }
public:
  prefixbuf(prefixer const& prefix, std::streambuf* sbuf)
      : sbuf(sbuf)
      , m_prefix(prefix)
      , need_prefix(true) {
  }
  virtual void setBuffer(std::streambuf* sb) {
    sbuf = sb;
  }
};

class oprefixstream
    : public virtual prefixbuf
    , public std::ostream
{
public:
  oprefixstream(prefixer const& prefix, std::ostream& out)
    : prefixbuf(prefix, out.rdbuf())
    , std::ios(static_cast<std::streambuf*>(this))
    , std::ostream(static_cast<std::streambuf*>(this)) {
  }
  void setStream(std::ostream& out) {
    sbuf = out.rdbuf();
  }
};

prefixer makePrefixer(std::string name) {
  return [name]() { std::stringstream ss;
    ss<<"["<<currentTime()<<"] "<<name<<": ";
    return ss.str();
  };
}

namespace LOG {
  using namespace std;

  static ofstream devnull("/dev/null");
  int verb(0);
  oprefixstream trace(makePrefixer("TRACE"), std::cout);
  oprefixstream debug(makePrefixer("DEBUG"), std::cout);
  oprefixstream info(makePrefixer("INFO"), std::cout);
  oprefixstream warn(makePrefixer("WARN"), std::cout);
  oprefixstream error(makePrefixer("ERROR"), std::cerr);
  oprefixstream fatal(makePrefixer("FATAL"), std::cerr);

  void setStream(int verbosity,
                 int threshold,
                 oprefixstream &stream,
                 std::ostream &output) {
    if (verbosity<threshold) {
      stream.setStream(devnull);
    } else {
      stream.setStream(output);
    }
  }
  
  std::string getVerbosity() {
    switch(verb) {
    case FATAL: return "FATAL";
    case ERROR: return "ERROR";
    case WARN: return "WARN";
    case INFO: return "INFO";
    case DEBUG: return "DEBUG";
    case TRACE: return "TRACE";
    default: return "QUIET";
    }
  }

  void setVerbosity(int verbosity) {
    verb = verbosity;
    setStream(verbosity,TRACE,LOG::trace,std::cout);
    setStream(verbosity,DEBUG,LOG::debug,std::cout);
    setStream(verbosity,INFO,LOG::info,std::cout);
    setStream(verbosity,WARN,LOG::warn,std::cout);
    setStream(verbosity,ERROR,LOG::error,std::cerr);
    setStream(verbosity,FATAL,LOG::fatal,std::cerr);
  }
}
