// when you lit up my cigarette, got out the taxi, asked me to sing that song that I wrote for you...

#include <iostream>
#include <string>

class Entity
{
private:
    std::string m_Name;
    int mutable m_DebugCount = 0;

public:
    const std::string &GetName() const
    {
        m_DebugCount++;
        return m_Name;
    }
};