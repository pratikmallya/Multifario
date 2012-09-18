/* 
 *  %W%
 *  %D% %T%
 * 
 *  author: Mike Henderson mhender@watson.ibm.com
 */

extern char *shOutputName;

#include <shInternal.h>
#include <string.h>
#include <stdlib.h>

#ifdef __cplusplus
 extern "C" {
#endif

void shSetOutputFormat(char *type)
 {
  shOutputFormat=(char*)realloc(shOutputFormat,(strlen(type)+1)*sizeof(char));
  strcpy(shOutputFormat,type);

  return;
 }

#ifdef __cplusplus
}
#endif
