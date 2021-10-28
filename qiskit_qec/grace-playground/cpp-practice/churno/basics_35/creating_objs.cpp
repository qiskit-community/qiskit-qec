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