/*===========================================================================*
   Cloned from AbcMessage.cpp
 *===========================================================================*/

#include "headers/BiqMessage.h"
#include <cstring>

typedef struct {
    BIQ_Message internalNumber;
    int externalNumber;              // or continuation
    char detail;
    const char * message;
} Biq_message;


// Once the message name is added to the enum BIQ_Message in the header file
// add it to this list by filling out the parameters of Biq_message
// still trying to figure out externalNumber and detail
// see BIQ_BOUND_END in BiqModel.cpp to see how to pass a message
// Curently BiqModel.cpp is the only place set up for messages
static Biq_message us_english[]=
{
    {BIQ_BOUND_END,1,1, "Bounding Complete:\n\tBound: %f BestVal: %f\n\tReason: %s\n\tItterations %d\n\tN Cuts: %d"},
    {BIQ_DUMMY_END, 999999, 0, ""}
};


/* Constructor */
BiqMessage::BiqMessage(Language language)
    :
    CoinMessages(sizeof(us_english) / sizeof(Biq_message))
{
    language_ = language;
    strcpy(source_, "Biq");
    Biq_message * message = us_english;

    while (message->internalNumber != BIQ_DUMMY_END) {
	CoinOneMessage oneMessage(message->externalNumber, message->detail,
				  message->message);
	addMessage(message->internalNumber, oneMessage);
	message++;
    }

    // now override any language ones

    switch (language) {

    default:
	message = NULL;
	break;
    }

    // replace if any found
    if (message) {
	while (message->internalNumber != BIQ_DUMMY_END) {
	    replaceMessage(message->internalNumber, message->message);
	    message++;
	}
    }
}
