#define LOGGING_LEVEL 3
#include "../src/include/logger.hpp"
#include <iostream>

int logger_test()
{    
    LOG_WARNING("Ich bin eine Warnung und eine echte ", 0, ".");
    LOG_INFO("Ich bin eine Falschinformation mit ", 2, " Zahlen darin!");
    LOG_DEBUG("Ich helfe beim Debugging in ", 2.2, " Prozent der ", "Faelle.");
    LOG_ERROR("Ich bin ein Fehler und dekonstruiere nur statische Objekte.");  
  
    // wird nie erreicht
    return 0;
}

int main()
{
    return logger_test();
}
