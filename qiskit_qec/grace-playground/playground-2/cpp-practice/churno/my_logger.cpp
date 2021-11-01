#include <iostream>

// use static as much as you can
class Log {
    // public variables
    public:
        const int LogLevelInfo = 0;
        const int LogLevelWarning = 1;
        const int LogLevelError = 2;

    // private var
    private:
        int m_LogLevel;

    // public methods
    public:
        void SetLevel(int level) {
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


        // virtual

        virtual std::string GetName() {
            return "Entity";
        }


};

class LitLog : public Log {
    private:
        std::string m_Name;

    public:

        LitLog(const std::string& name) : m_Name(name) {

        }

        std::string GetName() override {
            return  m_Name;
        }
};

int main() {

    LitLog log("Charnique");
    log.SetLevel(log.LogLevelWarning);
    log.Warn("You gently shut the door, then not a sound");
    std::cout << log.GetName() << std::endl;

}