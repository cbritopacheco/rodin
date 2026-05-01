// CoronaryArteryAlerts.h
#ifndef EXAMPLES_HEART_CORONARYARTERY_CORONARYARTERYALERTS_H
#define EXAMPLES_HEART_CORONARYARTERY_CORONARYARTERYALERTS_H

#include <Rodin/Alert.h>

namespace Rodin::Examples::Heart
{
  class ZeroDPrefix : public Alert::MessagePrefix<Alert::CyanT>
  {
    public:
      using Parent = Alert::MessagePrefix<Alert::CyanT>;

      ZeroDPrefix()
        : Parent("[0D]")
      {}
  };

  class ThreeDPrefix : public Alert::MessagePrefix<Alert::MagentaT>
  {
    public:
      using Parent = Alert::MessagePrefix<Alert::MagentaT>;

      ThreeDPrefix()
        : Parent("[3D]")
      {}
  };

  class SNESPrefix : public Alert::MessagePrefix<Alert::YellowT>
  {
    public:
      using Parent = Alert::MessagePrefix<Alert::YellowT>;

      SNESPrefix()
        : Parent("[SNES]")
      {}
  };

  class KSPPrefix : public Alert::MessagePrefix<Alert::YellowT>
  {
    public:
      using Parent = Alert::MessagePrefix<Alert::YellowT>;

      KSPPrefix()
        : Parent("[KSP]")
      {}
  };

  class TimingPrefix : public Alert::MessagePrefix<Alert::BrightGrayT>
  {
    public:
      using Parent = Alert::MessagePrefix<Alert::BrightGrayT>;

      TimingPrefix()
        : Parent("[Timing max/rank]")
      {}
  };

  using ZeroDInfo = Alert::PrefixedMessage<ZeroDPrefix>;
  using ThreeDInfo = Alert::PrefixedMessage<ThreeDPrefix>;
  using SNESInfo = Alert::PrefixedMessage<SNESPrefix>;
  using KSPInfo = Alert::PrefixedMessage<KSPPrefix>;
  using TimingInfo = Alert::PrefixedMessage<TimingPrefix>;
}

#endif
