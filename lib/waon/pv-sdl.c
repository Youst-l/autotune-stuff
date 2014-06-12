#include "SDL.h"
#include "SDL_mixer.h"

/* Mix_Music actually holds the music information.  */
Mix_Music *music = NULL;

int main(void) {

  SDL_Surface *screen;
  SDL_Event event;
  int done = 0;

  /* We're going to be requesting certain things from our audio
     device, so we set them up beforehand */
  int audio_rate = 22050;
  Uint16 audio_format = AUDIO_S16; /* 16-bit stereo */
  int audio_channels = 2;
  int audio_buffers = 4096;

  SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO);

  /* This is where we open up our audio device.  Mix_OpenAudio takes
     as its parameters the audio format we'd /like/ to have. */
  if(Mix_OpenAudio(audio_rate, audio_format, audio_channels, audio_buffers)) {
    printf("Unable to open audio!\n");
    exit(1);
  }

  /* If we actually care about what we got, we can ask here.  In this
     program we don't, but I'm showing the function call here anyway
     in case we'd want to know later. */
  Mix_QuerySpec(&audio_rate, &audio_format, &audio_channels);

  /* We're going to be using a window onscreen to register keypresses
     in.  We don't really care what it has in it, since we're not
     doing graphics, so we'll just throw something up there. */
  screen = SDL_SetVideoMode(320, 240, 0, 0);

  while(!done) {
    while(SDL_PollEvent(&event)) {
      switch(event.type) {
        case SDL_QUIT:
          done = 1;
          break;
        case SDL_KEYDOWN:
        case SDL_KEYUP:
          handleKey(event.key);
          break;
      }
    }

    /* So we don't hog the CPU */
    SDL_Delay(50);

  }

  /* This is the cleaning up part */
  Mix_CloseAudio();
  SDL_Quit();

}
