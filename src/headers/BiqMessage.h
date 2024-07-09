
#ifndef BiqMessage_H_
#define BiqMessage_H_

#include "CoinMessageHandler.hpp"

enum BIQ_Message
{
    BIQ_BOUND_HEAD,
    BIQ_BOUND_DATA,
    BIQ_DUMMY_END
};

class BiqMessage : public CoinMessages
{

public:
  BiqMessage(Language language=us_en);

};


#endif