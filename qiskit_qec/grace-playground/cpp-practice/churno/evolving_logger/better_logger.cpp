#include <iostream>

// use static as much as you can
class Log {


    public:

        // set up
        enum LogLevel : char {
          LogLevelError = 0, LogLevelWarning, LogLevelInfo
        };

        // init
        Log(LogLevel level) {
            m_LogLevel = level;

        }

        ~Log(){
            std::cout << "destroyed"<< std::endl;
        }

        // actionable methods
        void SetLevel(LogLevel level) {
            m_LogLevel = level;
        }

        void Info(const char* message) {
            if (m_LogLevel >= LogLevelInfo){
                std::cout << "[Info]]: " << message << std::endl;
            }
        }

        void Warn(const char* message) {
            if (m_LogLevel >= LogLevelWarning) {
                std::cout << "[WARNING]: " << message << std::endl;
            }

        }


        void Error(const char* message) {
            if (m_LogLevel >= LogLevelError){
                std::cout << "[ERROR]]: " << message << std::endl;
            }

        }

    // private var
    private:
        LogLevel m_LogLevel;


};



int main() {

    Log log(Log::LogLevelError);
    log.SetLevel(Log::LogLevelWarning);
    log.Warn("You gently shut the door, then not a sound");

    //std::cin.get();
}